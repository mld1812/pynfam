#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import numpy as np
import sys
import threading
from functools import partial
import datetime

def square_value(inp):
    return inp*inp

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
    with storage['lock']: # acquire lock of storage space
        # check storage again
        if not (storage.keys() >= {'length','results'} and isinstance(storage['length'], int) \
                and storage['length'] > 0 and isinstance(storage['results'], list)):
            print(f"Something wrong in master={master}")
            return
        # store result
        result = future.result()
        storage['results'].append(result)
        # send results back once all the results have been collected
        if len(storage['results']) >= storage['length']:
            print(storage['results'])
            comm.send(storage['results'], dest=master, tag=tag)
            storage['length'] = None
            storage['results'] = []

t_master2lead = 1000
t_lead2master = 2000
t_lead2worker = 3000
t_worker2lead = 4000
p_kill = -1000

comm_world = MPI.COMM_WORLD
rank = comm_world.Get_rank()
group = 0
if rank < 4:
    group = 1 #Split off first 4 as 'masters'

comm_master = comm_world.Split(group, 1)

new_comm_rank = comm_master.Get_rank()

print(f"Hello from process {rank}, in split comm: {new_comm_rank}")
"""
if rank < 4: #Masters
    #Send some tasks to be done to the lead worker (comm_world rank 4)
    task_list = np.arange(rank, rank + 40)
    comm_world.send(task_list, tag=t_master2lead, dest=4)
    #Listen for return message and then print it out.
    msg = comm_world.recv(source=4,tag=t_lead2master)
    print(f"Master {rank} received: {msg} at time: { datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")
    comm_master.Barrier() #Wait for all masters to receive return
    #send kill signal
    if new_comm_rank == 0:
        print("Sent kill")
        comm_world.send(p_kill, tag=t_master2lead, dest=4)

else: #Workers
    if __name__ == "__main__":
        with MPICommExecutor(comm_master, root=0) as executor:
            if executor is not None:
                status_ = MPI.Status()
                storages = {}
                while True: 
                    #listen for message from master
                    msg = comm_world.recv(source=MPI.ANY_SOURCE, tag=t_master2lead, status=status_)
                    if isinstance(msg, int) and msg == p_kill:
                        break
                    else:
                        master = status_.Get_source()
                        if master not in storages:
                            storages[master] = {'lock': threading.Lock()}
                        storage = storages[master]
                        with storage['lock']:
                            storage.update({'length': len(msg), 'results': []})
                            print(len(msg))
                        callback = partial(_processOneResult, master=master, tag=t_lead2master, storage=storage, \
                                       comm=comm_world)
                        for task in msg:
                            #Submit task to a worker
                            future = executor.submit(square_value, task)
                            future.add_done_callback(callback)
"""
"""
if __name__ == "__main__":
  #Try u sing MPICommExecutor.
  with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
    if executor is not None:
       future = executor.submit(abs, -42)
       assert future.result() == 42
       answer = set(executor.map(abs, [-42, 42]))
       assert answer == {42}

  size = MPI.COMM_WORLD.Get_size()
  rank = MPI.COMM_WORLD.Get_rank()
  if rank == 0:
    print(MPI.Query_thread())
  name = MPI.Get_processor_name()
  print_hello(rank, size, name)
  """