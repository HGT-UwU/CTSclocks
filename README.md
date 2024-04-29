---
title: "Introduction to CTSclocks"
author:
 - name: "Huige Tong"
   affiliation: "Shanghai Institute of Nutrition and Health, CAS"
   email: "tonghuige2021@sinh.ac.cn"
 - name: "Xiaolong Guo"
   affiliation: "Shanghai Institute of Nutrition and Health, CAS"
   email: "guoxiaolong2022@sinh.ac.cn"
 - name: "Andrew E Teschendorff"
   affiliation: "Shanghai Institute of Nutrition and Health, CAS"
   email: "andrew@sinh.ac.cn"
date: "2024-04-295"
package: CTSclocks
output:
  BiocStyle::html_document:
    toc_float: true
---

# Summary

CTSclocks is an R package which includes epigenetic clocks that are not confounded by cell-type heterogeneity and that can yield biological age estimates at cell-type resolution.

# Installation

To install (you will also need to install the dependency package presto):

```r
library(devtools)
devtools::install_github("HGT-UwU/CTSclocks")
```

# References

Huige Tong, Xiaolong Guo, Qi Luo and Andrew E Teschendorff. 2024. Cell-type specific epigenetic clocks to quantify biological age at cell-type resolution. Submitted.
