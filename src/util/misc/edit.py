#!/usr/bin/env python

from FileEditor import *

editor = FileEditor()
editor.setIsTest(False)

editor.setFilter(r'Util_Module')
editor.setOld(r'Util_Module')
editor.setNew(r'Misc_Module')
editor.editFiles(".", "*.h")
editor.editFiles(".", "*.mod")
