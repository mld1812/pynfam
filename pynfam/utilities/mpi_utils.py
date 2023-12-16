# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import str
# -------------- Utilities -----------------
import sys
import numpy as np
# ---------- MPI Implementation ------------
do_mpi = False
do_mpi_futures = False
try:
    from mpi4py import MPI
    if MPI.COMM_WORLD.Get_size()>1:
        do_mpi = True
        try:
            from mpi4py.futures import MPICommExecutor
            do_mpi_futures = True
        except ImportError:
            if MPI.COMM_WORLD.Get_rank()==0:
                print(u"*** mpi4py.futures not supported. ***")
    else:
        print(u"*** Only one MPI task requested. Performing serial calculation. ***")
except ImportError:
    print(u"*** mpi4py not supported. Performing serial calculation. ***")
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

# ------------------------------------------------------------------------------
# Tag codes (these are tags used for different data)
t_master2lead = 10 # master pinging lead
t_worker2lead = 20 # worker pinging lead
t_masterID    = 30 # lead sending master rank to worker
t_newtask     = 40 # master sending new task to worker
t_resulterr   = 50 # worker sending result error string (or initial ping) to master
t_resulttask  = 60 # worker sending result task to master

# Ping codes (these are sent/recieved)
p_kill =-9 # Terminate the program. This is compared with ranks, so use negative number
p_ping = 1 # Meaningless data sent as a ping. This is compared with a str, so any number is fine.

# ------------------------------------------------------------------------------
def _runOneTask(task, stdout=False, dbg=0):
    """
    Picklable wrapper function to return task and error.

    Args:
        task (fortProcess): A single task to run.
        stdout (bool): Print to stdout in real time if True.
        dbg (int): If !=0 don't actually run the exectuble.

    Returns:
        fortProcess: The input task, with now populated output.
        str, None: stderr string or None.
    """
    err_msg = task.runExe(stdout_bool=stdout, debug=dbg)
    return task, err_msg

# ------------------------------------------------------------------------------
def runtasks_master(task_list, comm_world, stdout=False, dbg=0):
    """
    Execute a list of fortran objects serially or assigning tasks to workers
    from a shared worker pool.

    Args:
        task_list (list of fortProcess): List of tasks.
        comm_world (MPI.Intracomm, int): The MPI communicator, or int if no MPI.
        stdout (bool): Print to stdout in real time if True.
        dbg (int): If !=0 don't actually run the exectuble.

    Returns:
        list of fortProcess: Input tasks with now populated output.
        bool: Was there an error?
        str, None: stderr string or None.
    """

    # Initialize
    finished, err_runs = [], []
    err_bool, err_msg = False, u''
    if not task_list:
        return finished, err_bool, err_msg


    # Allow serial runs by passing integer Comm
    try:
        rank, commsize = pynfam_mpi_traits(comm_world)
    except AttributeError:
        rank, commsize = 0, 1

    # Serial execution
    if commsize==1:
        for task in task_list:
            result_task, result_err = _runOneTask(task, stdout, dbg)
            finished.append(result_task)
            err_runs.append(result_err)

    # Dynamic parallel task allocation
    else:
        # Initialize
        unfinished_ids = list(range(len(task_list)))
        status_ = MPI.Status()
        task_sent = True

        # MASTER
        while True:
            # Ping lead worker requesting resources (NB!)
            # *** Should only send out as many pings as there are tasks.
            # *** Only send out a new ping once we've sent out the corresponding task. We might
            #     receive a result from another worker BEFORE we receive the response to this
            #     ping, in which case the task will never be sent and a worker will deadlock.
            if unfinished_ids and task_sent:
                comm_world.send(p_ping, dest=0, tag=t_master2lead)
                task_sent = False

            # Receive a ping or result from one of available workers
            result_err = comm_world.recv(source=MPI.ANY_SOURCE, status=status_, tag=t_resulterr)
            worker = status_.Get_source()

            # If we got a result, store it and loop, the worker is not looking for any further communication.
            if result_err != p_ping:
                result_task = comm_world.recv(source=worker, tag=t_resulttask)
                finished.append(result_task)
                err_runs.append(result_err)

            # If we got a worker ping and there are tasks, send one out
            elif unfinished_ids:
                task_id = unfinished_ids.pop()
                comm_world.send(task_list[task_id], dest=worker, tag=t_newtask)
                task_sent = True

            if len(finished)>=len(task_list):
                break

    errors = [e for e in err_runs if e is not None]
    if errors:
        err_bool = True
        err_msg = errors[0]

    return finished, err_bool, err_msg

# ------------------------------------------------------------------------------
def runtasks_worker(comm_world, comm_worker, stdout=False, dbg=0):
    """
    The worker pool in charge of running fortProcess objects for load balancing
    MPI calculations. A lead worker assigns resources, while workers recv/send
    fortProcess tasks/results.

    Args:
        comm_world (MPI.Intracomm, int): The MPI communicator COMM_WORLD
        comm_worker (MPI.Intracomm, int): The MPI communicator for workers
        dbg (int): If !=0 don't actually run the exectuble.

    Returns:
        list of fortProcess: Input tasks with now populated output.
        bool: Was there an error?
        str, None: stderr string or None.
    """

    # Initialize variables
    rank_world = comm_world.Get_rank()
    rank_worker = comm_worker.Get_rank()
    nr_workers = comm_worker.Get_size() - 1
    kill_total = 0
    status_ = MPI.Status()

    # LEAD WORKER (rank_worker=0 is rank_world=0)
    if rank_world == 0:
        while True:
            # Look for master requesting resources
            mping = comm_world.recv(source=MPI.ANY_SOURCE, status=status_, tag=t_master2lead)
            master = status_.Get_source()
            # Look for available worker
            wping = comm_worker.recv(source=MPI.ANY_SOURCE, status=status_, tag=t_worker2lead)
            worker = status_.Get_source()
            # Send the master rank to the worker (or signal to terminate)
            if mping != p_kill:
                comm_worker.send(master, dest=worker, tag=t_masterID)
            else:
                comm_worker.send(p_kill, dest=worker, tag=t_masterID)
                kill_total += 1
            if kill_total == nr_workers:
                break
    # WORKER POOL
    else:
        while True:
            # Ping the leader to recieve rank of master in need of resources
            comm_worker.send(p_ping, dest=0, tag=t_worker2lead)
            master = comm_worker.recv(source=0, tag=t_masterID)

            # Ping the master, get task, run it, send back result
            if master != p_kill:
                comm_world.send(p_ping, dest=master, tag=t_resulterr)
                task = comm_world.recv(source=master, tag=t_newtask)

                result_task, result_err = _runOneTask(task, stdout)

                comm_world.send(result_err, dest=master, tag=t_resulterr)
                comm_world.send(result_task, dest=master, tag=t_resulttask)
            else:
                break

# ------------------------------------------------------------------------------
def runtasks_killsignal(comm_world, comm_master):
    '''
    When the program is finished, send a signal to all workers in the worker
    pool to break their loop and exit.

    Args:
        comm_world (MPI.Intracomm, int): The MPI communicator COMM_WORLD
        comm_world (MPI.Intracomm, int): The MPI communicator for masters
    '''
    # Sync up masters before beginning the finalize process
    comm_master.Barrier()

    # nr_workers = total_processes - nr_masters - 1_lead_worker
    rank_master, commsize_master = pynfam_mpi_traits(comm_master)
    commsize_world = comm_world.Get_size()
    nr_workers = commsize_world - commsize_master - 1

    # Send 1 kill signal for each worker to the lead
    if rank_master == 0:
        for i in range(nr_workers):
            comm_world.send(p_kill, dest=0, tag=t_master2lead)

# ------------------------------------------------------------------------------
def pynfam_mpi_traits(Comm):
    """
    Get rank and size of a communicator. Dummy communicator for serial
    calculations returns rank 0 size 1.

    Args:
        Comm (MPI.Intracomm): The communicator.

    Returns:
        int: Rank.
        int: Size.
    """
    if do_mpi:
        rank     = Comm.Get_rank()
        commsize = Comm.Get_size()
    else:
        rank     = 0
        commsize = 1
    return rank, commsize

# ------------------------------------------------------------------------------
def pynfam_comm_split(Comm, nr_calcs, min_mpi_tasks=1, key=1):
    """
    Split a communicator into groups based on total number of calculations
    and minimum number of tasks in any group.

    Args:
        Comm (MPI.Intracomm): The communicator.
        nr_calcs (int): Number of calculations to be performed.
        mpi_mpi_tasks (int): Minimum number of tasks in any group (default 1).
        key (int): Order ranks in newcomm, tie retains order from original comm. (default 1).

    Returns:
        int: Number of groups.
        int: Group label.
        MPI.Intracomm: The new communicator.
    """
    if do_mpi:
        rank, comm_size = pynfam_mpi_traits(Comm)
        # Integer division here
        nr_groups = nr_calcs
        nr_even_split = comm_size//nr_groups
        # If nr_even_split = 0, min_mpi_tasks>=1 fixes this
        if min_mpi_tasks > nr_even_split:
            nr_groups = comm_size//min_mpi_tasks
        # Split the communicator
        group   = rank%nr_groups
        Newcomm = Comm.Split(group, key)
    else:
        nr_groups, group, Newcomm = 1, 0, 0

    return nr_groups, group, Newcomm

# ------------------------------------------------------------------------------
def pynfam_gather(dist_lists, Comm, root):
    """
    Gather and flatten distributed lists of results to rank=root.

    Args:
        dist_lists (list): List of data to be gathered.
        Comm (MPI.Intracomm): The communicator.
        root (int): The rank recieving the data.

    Returns:
        list: Gathered data on rank=root.
    """
    if do_mpi:
        rank = pynfam_mpi_traits(Comm)[0]
        gath_list, fin_list = [], []
        gath_list = Comm.gather(dist_lists, root=root) # list of lists
        # flatten, skipping empties
        if rank == root:
            fin_list = [item for sublist in gath_list for item in sublist]
    else:
        fin_list = dist_lists # only 1 list if no mpi

    return fin_list

# ------------------------------------------------------------------------------
def pynfam_abort(Comm=None, msg=u"Abort was called."):
    """
    MPI compatible program exit with message. Abort kills all teams/tasks.

    Args:
        Comm (MPI.Intracomm): The communicator (default None).
        msg (str): Message to display before exiting the program.
    """
    if not isinstance(msg,list): msg = [msg]

    print()
    print(u" ************************************************** ")
    print(u"  ERROR raised:")
    for m in msg: print(u"    "+str(m))
    print(u" ************************************************** ")
    print()
    sys.stdout.flush()

    if do_mpi:
        Comm.Abort()
    else:
        exit()

# ------------------------------------------------------------------------------
def pynfam_warn(msg, group=u'', error=None, Comm=None):
    """
    MPI compatible non-fatal warning message. Checks for any error on any task.
    If error is detected, print message and broadcast error boolean to all tasks.
    If only msg provided, rank=0 prints message and returns True (None on others).

    Args:
        msg (str): Warning message to display.
        group (str): Group label on which the warning is raised (default '')
        error (bool): Error indicator (default None)
        Comm (MPI.Intracomm, int): The communicator, or int if no mpi.

    Returns:
        bool
    """
    # If error is None, default error = True on all tasks
    if error is None: error=True

    # Gather the errors at rank 0. This is a blocking communication.
    if do_mpi and Comm is not None:
        gerr_list = None
        gerr_list = Comm.gather(error, root=0)
        rank = pynfam_mpi_traits(Comm)[0]
    else:
        gerr_list = [error]
        rank = 0

    # If any errors are true, skip = true and print message
    if not isinstance(msg,list): msg = [msg]
    skip = None
    if rank == 0:
        skip = any(e for e in gerr_list)
        if skip:
            print()
            print(u" ************************************************** ")
            print(u"  WARNING raised:")
            if group: print(u"    Warning in pynfam calc "+str(group))
            for m in msg: print(u"    "+str(m))
            print(u" ************************************************** ")
            print()
            sys.stdout.flush()
    if do_mpi and Comm is not None: skip = Comm.bcast(skip, root=0)

    return skip

# ------------------------------------------------------------------------------
sys_excepthook = sys.excepthook
def mpi_excepthook(v, t, tb):
    """
    Extend sys.excepthook to call mpi abort in the case of unhandled
    excpetions being raised. This prevents deadlock. Also print stderr.
    """
    print(file=sys.stderr)
    print(u" ************************************************** ",file=sys.stderr)
    print(u"  Unhandled Exception Raised:",file=sys.stderr)
    sys_excepthook(v, t, tb)
    print(u" ************************************************** ",file=sys.stderr)
    print(file=sys.stderr)
    sys.stderr.flush()
    if do_mpi: MPI.COMM_WORLD.Abort(1)
sys.excepthook = mpi_excepthook
