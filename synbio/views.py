from . import view

class GenericView(view.BaseView):
    def _validate(self):
        pass # allow all clumps

    @view.view_method
    def count_sequences(self):
        return len(self._clump._sequences)
