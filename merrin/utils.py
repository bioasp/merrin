
import time
import threading

def count_and_time(func):
    def helper(self, *args, **kwargs):
        attr = f"__stats_{func.__name__}"
        if not hasattr(self, attr):
            stats = {"calls": 0, "total_duration": 0}
            setattr(self, attr, stats)
        else:
            stats = getattr(self, attr)
        stats["calls"] += 1
        def get_time():
            #return time.time()
            thread_id = threading.get_ident()
            clk = time.pthread_getcpuclockid(thread_id)
            return time.clock_gettime(clk)
        t0 = get_time()
        res = func(self, *args, **kwargs)
        stats["total_duration"] += get_time() - t0
        return res
    helper.calls = 0
    helper.total_duration = 0
    helper.__name__= func.__name__
    return helper


