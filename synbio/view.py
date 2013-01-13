class ViewError(Exception):
    pass

class BaseView:
    def __init__(self, clump):
        self._clump = clump
        self._validate()

def view_method(func):
    def out(self_view, *args, **kwargs):
        self_view._validate()
        return func(self_view, *args, **kwargs)

    return out
