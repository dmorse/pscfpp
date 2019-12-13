#!/usr/bin/env python

from pscfpp.FileEditor import *

editor = FileEditor()
editor.setIsTest(False)

editor.setFilter(r'wFieldGrids')
editor.setOld(r'wFieldGrids')
editor.setNew(r'wFieldsRGrid')
editor.editFiles("pspc", "*")
 
editor.setFilter(r'wFieldGrid')
editor.setOld(r'wFieldGrid')
editor.setNew(r'wFieldRGrid')
editor.editFiles("pspc", "*")
 
editor.setFilter(r'cFieldGrids')
editor.setOld(r'cFieldGrids')
editor.setNew(r'cFieldsRGrid')
editor.editFiles("pspc", "*")

editor.setFilter(r'cFieldGrid')
editor.setOld(r'cFieldGrid')
editor.setNew(r'cFieldRGrid')
editor.editFiles("pspc", "*")
 
editor.setFilter(r'wFieldDfts')
editor.setOld(r'wFieldDfts')
editor.setNew(r'wFieldsKGrid')
editor.editFiles("pspc", "*")
 
editor.setFilter(r'wFieldDft')
editor.setOld(r'wFieldDft')
editor.setNew(r'wFieldKGrid')
editor.editFiles("pspc", "*")
 
editor.setFilter(r'cFieldDfts')
editor.setOld(r'cFieldDfts')
editor.setNew(r'cFieldsKGrid')
editor.editFiles("pspc", "*")

editor.setFilter(r'cFieldDft')
editor.setOld(r'cFieldDft')
editor.setNew(r'cFieldKGrid')
editor.editFiles("pspc", "*")
 
