# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

from ..Util.Observer import Observable


class Vector(Observable):
    _Args = {"x": float("nan"), "y": float("nan"), "z": float("nan")}
    _FORMAT_STRING = "({0.x},{0.y},{0.z})"

    def __init__(self, *args):
        if len(args) == 0:
            nan = float("nan")
            self.x = nan
            self.y = nan
            self.z = nan
        elif len(args) == 1:
            self.x, self.y, self.z = args[0]
        elif len(args) == 3:
            self.x, self.y, self.z = args
        else:
            raise ValueError

        return

    def __repr__(self):
        return "Vector" + str(self)

    def __str__(self):
        return self._FORMAT_STRING.format(self)
