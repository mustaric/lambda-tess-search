# -*- coding: utf-8 -*-

"""
Created on Mon Aug  3 13:58:38 2020

TODO
x I need a timeout for each task
x I need an error_handler to return None if a task fails
o I need more control over maxtask, chunksize, etc.
@author: fergal
"""



import concurrent.futures
import multiprocessing
import functools
import itertools
import traceback
import datetime
import asyncio
import shutil
import sys


def default_error_response(func, task, exc):
    """Helper function"""
    return None



def parmap(func, *args, fargs=None,
               engine='multi',
               timeout_sec=1e6,
               n_simul=None,
               on_error=default_error_response,
            ):
    """Apply map in parallel.

    This is a python equivalent to Matlab's parfor function. Runs
    trivially parallisable problems in multiple parallel processes

    Examples
    ----------

    .. code-block::

        x = np.arange(5)

        #Parallise a call to a function with one argument
        def sqr(x):
            return x*x
        parmap(sqr, x)
        >>> [0, 1, 4, 9, 16, 25]

        #Parallelise a function with two arguments
        def hypot(x, y):
            return np.sqrt(x**2 + y**2)

        #hypot is called on every combination of x[i] and x[j].
        #result has twenty-five elements
        parmap(hypot, x, x)
        >>> [0, 1, 2, ... 7.071]

        #Parallelise a function with  a configuration option
        def power(x, n):
            return x**n

        #parmap works accepts both positional and keyword arguments as keyword arguments
        function_args = dict(n=3)
        result = parmap(hypot, x, fargs=function_args)

        def hypotn(x, y, n):
            return x**n + y**n
        result = parmap(hypot, x, x, fargs=function_args)


    Arguments
    -----------
    func
        The function to apply. Decorated functions may raise errors
        due to the way Python passes the function to the child processes

    *args
        One or more iterables. `func` is applied to every combination
        of elements. For example, if args = [(1,2), (3,4)], the func
        is called 4 times, with arguments of (1,3), (1,4), (2,3), (2,4)
    fargs
        (dict) A dictionary of non-iterable arguments. See examples above
    timeout_sec
        (int) Kill a task if it does not return after this many seconds.
        The default value is 11 days. If your tasks take that long you
        need to refactor your code.
    on_error
        (function) Error handling function. See below.
    n_simul
        (int) Number of simultaneous tasks to process.
    engine
        (str) Select an engine to use for your processing.

        multi
            (Default) Use the multiprocessing module
        serial
            Run in a single thread on a single process. This is a useful
            debugging mode.
        threads
            Use parallel threads in a single process. Due to python's GIL,
            only one thread can execute at a time. However, threads
            use less memory than multiprocessing. Threads may be more
            suitable if your code spends a lot of time waiting on IO.
        async
            Experimental. Use asyncio's concurrent processing.
            May be lower memory than threads, but `pfunc` must
            be an async function

        'multi' is a good default. Threads and async may give
        your lower memory requirements for IO heavy tasks. If you
        select the 'async' engine, `func` must be an async function.

    Returns
    ---------
    A list of return values from the function.

    Notes
    ------------------
    Error Handling
        Error handling depends on the engine used. serial and async throw
        an exception and let you deal with it. This is invaluable for
        debugging.

        Debugging errors in parallel processes is more difficult. For
        'multi' and 'threads' engines, parmap replaces the result for any failed
        task with the result of calling the function specified by `on_error`.
        The default error handler returns **None**, but you can write your own error
        handler to return what you like.

        Bear in mind that the error handler is not called in a separate processor,
        and is expected to return quickly. See `default_error_response()` for
        the signature of the error handler

    Small Tasks
        For every short tasks (such as those shown in the example) the cost
        of forking a new process can exceed the runtime of the process. The
        'multi' engine only results in speed gains for functions that take
        more than about half a second to run.

    Lambda Functions
        Lambda functions can't be passed to the 'multi' or 'thread' engines
        because of internal details of how python tries (or fails) to pickle
        these functions. Internal functions (i.e functions defined inside
        other functions) may also cause trouble, and aren't well tested.

    """
    if fargs is None:
        fargs = {}

    tasks = list(itertools.product(*args))
    if engine != 'serial':
        #Not needed for the serial engine
        func = BacktraceCatcher(func)
    pfunc = functools.partial(func, **fargs)

    engine_dict = load_engine_dict()
    try:
        engine = engine_dict[engine]
    except KeyError:
        raise KeyError("Unrecognised engine %s. Must be one of %s" %(engine, engine_dict.keys()))

    pfunc.__name__ = func.__name__
    reporter='bar'
    reporter = get_reporter(reporter)  #Progress reporter
    results = engine(pfunc, tasks, n_simul, timeout_sec, on_error, reporter)
    return results

def get_reporter(request):

    opts = dict(
        silent=NullReporter(),
        text=TextReporter(),
        bar=ProgressBarReporter(),
    )

    if isinstance(request, str):
        try:
            return opts[request]
        except KeyError:
            return KeyError("Requested reporter must be one of %s" %(opts.keys()))

    #Assume request is a Reporter like object
    return request        
class BacktraceCatcher():
    """Stuff the backtrace into the exception value

    Backtrace objects can't be pickled, so a multi-processing task
    that fails can't propegate the backtrace back to the main process.

    This wrapper converts the backtrace to a string and inserts it into
    the exception to work around this limitation.

    This helps us by making the exceptions created during multiprocessing
    easier to debug.
    """
    def __init__(self, func):
        self.func = func
        self.__name__ = func.__name__

    def __call__(self, *args, **kwargs):
        try:
            return self.func(*args, **kwargs)
        except Exception as exc:
            etype, evalue, etrace = sys.exc_info()
            msg = traceback.format_tb(etrace)
            exc.args += tuple(msg)
            raise exc


def load_engine_dict():
    engine_dict = { 
        'multi':parallel_apply,
        'threads':thread_apply,
        'async':async_apply,
        'serial': linear_apply,
    }
    return engine_dict


def linear_apply(pfunc, tasks, n_simul, timout_sec, on_error, reporter):
    """The simplest engine. Run tasks serially"""
    results = []
    reporter(-1, len(tasks))  #Shows 0 tasks complete
    for i, task in enumerate(tasks):
        results.append(pfunc(*task))
        reporter(i, len(tasks))
    return results


def thread_apply(pfunc, tasks, n_thread, timeout_sec, on_error, reporter):
    """Run each task in a thread"""
    if n_thread is None:
        n_thread = 5

    results = []
    reporter(-1, len(tasks))  #Shows 0 tasks complete
    with concurrent.futures.ThreadPoolExecutor(n_thread) as executor:
        future_list = map(lambda x: executor.submit(pfunc, *x), tasks)
        future_list = list(future_list)

        for i in range(len(future_list)):
            fut = future_list[i]
            try:
                r = fut.result(timeout=timeout_sec)
            except Exception as exc:
                warn_on_error(pfunc, tasks[i], exc)
                r = on_error(pfunc, tasks[i], exc)

            reporter(i, len(future_list))
            results.append(r)
    return results


def parallel_apply(pfunc, tasks, n_cpu, timeout_sec, on_error, reporter):
    """Run each task in a separate process"""
    if n_cpu is None:
        n_cpu = multiprocessing.cpu_count() - 1

    results = []
    reporter(-1, len(tasks))  #Shows 0 tasks complete
    with multiprocessing.Pool(n_cpu, maxtasksperchild=1) as pool:
        process_list  = map(lambda x: pool.apply_async(pfunc, x), tasks)
        process_list = list(process_list)  #Start the tasks running

        for i in range(len(process_list)):
            p = process_list[i]
            try:
                r = p.get(timeout=timeout_sec)
            except Exception as exc:
                warn_on_error(pfunc, tasks[i], exc)
                r = on_error(pfunc, tasks[i], exc)

            reporter(i, len(process_list))
            results.append(r)
    return results


def async_apply(pfunc, tasks, n_simul, timeout_sec, on_error, reporter):
    return asyncio.run(_async_for(pfunc, tasks, timeout_sec, reporter))

async def _async_for(pfunc, tasklist, timeout_sec, reporter):
    futurelist = []
    for task in tasklist:
        future = asyncio.create_task(pfunc(*task))
        futurelist.append(future)

    out = []
    for f in futurelist:
        res = await f
        out.append(res)
    return out


def warn_on_error(func, task, error):
    msg = "WARN: Function %s on task %s\n failed with error: '%s'"
    msg = msg % (func.__name__, str(task), repr(error))
    print(msg)




#testing code
import numpy as np
import pytest

def sqr(x):
    print("The squre of %i is %i" %(x, x**2))
    return x*x


def power(x, n):
    return x**n


def hypot(x, y):
    return np.sqrt(x**2 + y**2)

def hypotn(x, y, n):
    return x**n + y**n


async def asqr(x):
    await asyncio.sleep(2)
    return x**2


def test_sqr():
    x = np.arange(10)

    for sp in [True, False]:
        res = parmap(sqr, x, single_process=sp)
        assert np.all(res == x **2)


def test_pow():
    x = np.arange(10)

    for sp in [False ]: #, False]:
        res = parmap(power, x, fargs={'n':3}, single_process=sp)
        assert np.all(res == x **3)

def test_hypot():
    x = np.arange(10)
    y = 10 + np.arange(10)

    for sp in [True, False]:
        res = parmap(hypot, x, y, single_process=sp)
        assert len(res) == 100
    return res


def test_hypotn():
    x = np.arange(10)
    y = np.arange(10)

    for sp in [True, False]:
        res = parmap(hypotn, x, y, fargs={'n':3}, single_process=sp)
        assert len(res) == 100
        assert res[-1] == 2*9**3


def failing_task(x):
    if x == 2:
        1/0
    return x

async def afailing_task(x):
    if x == 2:
        1/0
    return x


def test_task_that_fails():
    x = np.arange(5)

    with pytest.raises(ZeroDivisionError):
        res = parmap(failing_task, x, single_process=True)

    res = parmap(failing_task, x, single_process=False)
    assert res[2] is None
    assert res[1] == 1


def test_async():
    x = np.arange(10)
    res = parmap(asqr, x, engine='async')
    assert np.all(res == x**2)

def test_async_fail():
    x = np.arange(10)

    with pytest.raises(ZeroDivisionError):
        parmap(afailing_task, x, engine='async')


class NullReporter():
    """Don't report anything on iteration"""
    def __call__(self, itr, num):
        pass 


class TextReporter():
    def __init__(self):
        self.start_time = datetime.datetime.now()

    def __call__(self, itr, num):
        now = datetime.datetime.now()
        elapsed = now - self.start_time 
        elapsed_sec = elapsed.total_seconds()
        print("%i/%i tasks complete (%i sec elapsed)" %(itr+1, num, elapsed_sec))



class ProgressBarReporter():
    """A Tqdm-like progress bar.
    
    I tried using tqdm, but this usecase is just a bit
    outside it's comfortzone
    """
    def __init__(self):
        self.start_time = datetime.datetime.now()

    def __call__(self, itr, num):
        itr += 1

        # bar = self.create_bar(itr, num)
        elapsed = self.compute_elapsed(itr, num)
        msg = "%i/%i %s " %(itr, num, elapsed)

        bar = self.create_bar(itr, num, len(msg))
        msg = bar + msg
        print(msg)

    def create_bar(self, itr, num, nchars):
        width = shutil.get_terminal_size().columns - nchars - 10

        frac = int(width * itr/num)
        remain = width - frac
        #Not sure which unicode character looks better
        # bar = "[" + "█" * frac + "·" * remain + "] "
        # bar = "[" + "■" * frac + "·" * remain + "] "
        
        bar = "[" + "▬" * (frac-1) + "▶" + "·" * remain + "] "
        
         	 
        return bar
    
    def compute_elapsed(self, itr, num):
        now = datetime.datetime.now()
        elapsed = now - self.start_time 
        elapsed_sec = elapsed.total_seconds()
        return "(%i sec)" %(elapsed_sec)

