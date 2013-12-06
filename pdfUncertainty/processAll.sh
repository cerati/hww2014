#!/bin/bash

# low mass

root -b -q process2D.C+\(\"analysis22_125\",\"qqWW_DF_0j\",\"cteq6ll_0\"\)
root -b -q process2D.C+\(\"analysis22_125\",\"qqWW_DF_1j\",\"cteq6ll_0\"\)
root -b -q process2D.C+\(\"analysis22_125\",\"qqWW_DF_VBF\",\"cteq6ll_0\"\)

root -b -q process2D.C+\(\"analysis22_125\",\"ggWW_DF_0j\",\"cteq6ll_0\"\)
root -b -q process2D.C+\(\"analysis22_125\",\"ggWW_DF_1j\",\"cteq6ll_0\"\)

# high mass

root -b -q process2D.C+\(\"analysis21_500\",\"qqWW_DF_0j\",\"cteq6ll_0\"\)
root -b -q process2D.C+\(\"analysis21_500\",\"qqWW_DF_1j\",\"cteq6ll_0\"\)

root -b -q process2D.C+\(\"analysis21_500\",\"ggWW_DF_0j\",\"cteq6ll_0\"\)
root -b -q process2D.C+\(\"analysis21_500\",\"ggWW_DF_1j\",\"cteq6ll_0\"\)

# combine...

hadd results2/PDFUncertainty_LowMass.root \
    results2/PDFUncertainty_analysis22_125_ggWW_DF_0j.root results2/PDFUncertainty_analysis22_125_ggWW_DF_1j.root \
    results2/PDFUncertainty_analysis22_125_qqWW_DF_0j.root results2/PDFUncertainty_analysis22_125_qqWW_DF_1j.root 

hadd results2/PDFUncertainty_HighMass.root \
    results2/PDFUncertainty_analysis21_500_ggWW_DF_0j.root results2/PDFUncertainty_analysis21_500_ggWW_DF_1j.root \
    results2/PDFUncertainty_analysis21_500_qqWW_DF_0j.root results2/PDFUncertainty_analysis21_500_qqWW_DF_1j.root   

