# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import str
# -------------- Utilities -----------------
import sys
import time
import numpy as np
from functools import partial
import threading
# ---------- MPI Implementation ------------
do_mpi = False
try:
    from mpi4py import MPI
    from mpi4py.futures import MPICommExecutor
    do_mpi = True
except ImportError:
    print(u"*** mpi4py not supported. Performing serial calculation. ***")
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

# ------------------------------------------------------------------------------
# Tag code
t_master2lead = 10000 # master sending tasks to lead

# Ping code
p_kill =-9 # Terminate the program.

# Type of thread lock
type_lock = type(threading.Lock())

# ------------------------------------------------------------------------------
def _comm_recv(comm, source, tag, status=None):
    """
    This function replaces comm.recv() to avoid 100% CPU consumption.
    The implementation is based on following code:
    https://github.com/mpi4py/mpi4py/blob/master/src/mpi4py/futures/_core.py

    Args:
        comm (MPI.Intracomm): The MPI communicator.
        source (int): The rank of source process where the message is sent.
        tag (int): The tag of the MPI communication.
        status (MPI.Status, None): The status of the MPI communication.

    Returns:
        None
    """
    tval = 0.0; tmax = 0.001; tmin = 0.001 / (1 << 10)
    if not isinstance(status, MPI.Status):
        status = MPI.Status()
    while not comm.iprobe(source=source, tag=tag, status=status):
        time.sleep(tval)
        tval = min(tmax, max(tmin, tval*2))
    source, tag = status.source, status.tag
    message = comm.recv(source=source, tag=tag, status=status)
    return message

# ------------------------------------------------------------------------------
def _runOneTask(task, stdout=False, dbg=0):
    """
    Picklable wrapper function to return task and error.

    Args:
        task (fortProcess): A single task to run.
        stdout (bool): Print to stdout in real time if True.
        dbg (int): If !=0 don't actually run the executable.

    Returns:
        fortProcess: The input task, with now populated output.
        str, None: stderr string or None.
    """
    try:
        err_msg = task.runExe(stdout_bool=stdout, debug=dbg)
        return task, err_msg
    except BaseException as ex:
        return task, str(ex)

# ------------------------------------------------------------------------------
def _processOneResult(future, master, tag, storage, comm):
    """
    Callback function called by lead worker to store result returned by a non-lead
    worker, and send results back to the corresponding master process once all the 
    results have been collected. 

    Args:
        future (concurrent.futures.Future): The Future obj of the task.
        master (int): The rank of master process to send result back.
        tag (int): The tag to use in MPI communication.
        storage (dict): The storage space with following components:
            storage['length'] (int): The total number of results to process.
            storage['results'] (list): Space to store results.
            storage['lock'] (threading.Lock): Thread lock to avoid simultaneous
                                              modifications on storage. 
        comm (MPI.Intracomm): The MPI communicator.

    Returns:
        None.
    """
    error_msg = 'Storage space for index={:} is not properly initialized; MPI will hang.'
    # check storage
    if not (isinstance(storage, dict) and 'lock' in storage and isinstance(storage['lock'], type_lock)):
        pynfam_warn(error_msg.format(tag)); return
    with storage['lock']: # acquire lock of storage space
        # check storage again
        if not (storage.keys() >= {'length','results'} and isinstance(storage['length'], int) \
                and storage['length'] > 0 and isinstance(storage['results'], list)):
            pynfam_warn(error_msg.format(tag)); return
        # store result
        result = future.result()
        storage['results'].append(result)
        # send results back once all the results have been collected
        if len(storage['results']) >= storage['length']:
            comm.send(tuple(zip(*storage['results'])), dest=master, tag=tag)
            storage['length'] = None
            storage['results'] = []

# ------------------------------------------------------------------------------
def runtasks_master(task_list, comm_admin, stdout, index, dbg=0):
    """
    Execute a list of fortran objects serially or assigning tasks to workers
    from a shared worker pool.

    Args:
        task_list (list of fortProcess): List of tasks.
        comm_admin (MPI.Intracomm, int): The MPI communicator for administrators 
                                         (lead worker and masters), or int if no MPI.
        stdout (bool): Print to stdout in real time if True.
        index (int): Index of the calc, used for tag in MPI communication.
        dbg (int): If !=0 don't actually run the executable.

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
        rank, commsize = pynfam_mpi_traits(comm_admin)
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
        # MASTER sends index and tasks to LEAD WORKER
        comm_admin.send((index, task_list), dest=0, tag=t_master2lead)

        # MASTER waits for results from LEAD WORKER
        finished, err_runs = _comm_recv(comm_admin, source=0, tag=index)

    errors = [e for e in err_runs if e is not None]
    if errors:
        err_bool = True
        err_msg = errors[0]

    return list(finished), err_bool, err_msg

# ------------------------------------------------------------------------------
def runtasks_worker(comm_admin, comm_worker, stdout, dbg=0):
    """
    The worker pool in charge of running fortProcess objects for load balancing
    MPI calculations. A lead worker assigns resources, while workers recv/send
    fortProcess tasks/results.

    Args:
        comm_admin (MPI.Intracomm, int): The MPI communicator for administrators
                                         (lead worker and masters)
        comm_worker (MPI.Intracomm, int): The MPI communicator for workers
        stdout (bool): Print to stdout in real time if True.
        dbg (int): If !=0 don't actually run the executable.

    Returns:
        None.
    """

    # WORKER POOL
    with MPICommExecutor(comm_worker, root=0) as executor:
        if executor is not None:
            status_ = MPI.Status()
            storages = {}
            # LEAD WORKER receives tasks from masters and submits them to the executor
            while True:
                msg = _comm_recv(comm_admin, source=MPI.ANY_SOURCE, tag=t_master2lead, status=status_)
                if msg != p_kill: # Submit tasks
                    # Get index, tasks, master (source of message)
                    index, tasks = msg
                    master = status_.Get_source()
                    # Prepare storage space for results
                    if index not in storages:
                        storages[index] = {'lock': threading.Lock()}
                    storage = storages[index]
                    with storage['lock']:
                        storage.update({'length': len(tasks), 'results': []})
                    # Callback to store result and send results back to master
                    callback = partial(_processOneResult, master=master, tag=index, storage=storage, \
                                       comm=comm_admin)
                    for task in tasks: # Submit and add callback
                        future = executor.submit(_runOneTask, task, stdout, dbg)
                        future.add_done_callback(callback)
                else: # Wait all workers and shutdown executor
                    break

# ------------------------------------------------------------------------------
def runtasks_killsignal(comm_admin, comm_master):
    '''
    When the program is finished, send a signal to all workers in the worker
    pool to break their loop and exit.

    Args:
        comm_admin (MPI.Intracomm, int): The MPI communicator for administrators
                                         (lead worker and masters)
        comm_master (MPI.Intracomm, int): The MPI communicator for masters
    '''
    # Sync up masters before beginning the finalize process
    comm_master.Barrier()

    rank_master, commsize_master = pynfam_mpi_traits(comm_master)

    # Send kill signal
    if rank_master == 0:
        comm_admin.send(p_kill, dest=0, tag=t_master2lead)

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

    print(file=sys.stderr)
    print(u" ************************************************** ", file=sys.stderr)
    print(u"  ERROR raised:", file=sys.stderr)
    for m in msg: print(u"    "+str(m), file=sys.stderr)
    print(u" ************************************************** ", file=sys.stderr)
    print(file=sys.stderr)
    sys.stderr.flush()

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