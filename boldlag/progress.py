"""Progress reporting hook.

Library code calls :func:`report(stage, fraction)`; front ends register a callback
with :func:`set_callback`.  Without a callback nothing happens.
"""
_cb = None


def set_callback(cb):
    """``cb(stage: str, fraction: float)`` with fraction in 0..1 (overall job progress)."""
    global _cb
    _cb = cb


def report(stage, fraction):
    if _cb is not None:
        try:
            _cb(stage, max(0.0, min(1.0, float(fraction))))
        except Exception:       # a broken front end must not kill the job; BaseException (e.g. Streamlit's Stop) propagates
            pass


class Span:
    """Map sub-progress 0..1 of one stage onto [start, stop] of the overall job."""
    def __init__(self, stage, start, stop):
        self.stage, self.start, self.stop = stage, start, stop
    def __call__(self, frac):
        report(self.stage, self.start + (self.stop - self.start) * frac)
