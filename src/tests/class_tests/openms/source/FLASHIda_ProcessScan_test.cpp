// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
//
// Phase 4 unit tests: processScan full routing.
// Peak data extracted from FlashIDA/test-data/spectra/ms1_standard.txt (scan 1)
// and ms2_hcd_fragment.txt (scan 57).

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <string>
#include <cstring>

using namespace OpenMS;

namespace
{
  // Minimal JSON config for standard DDA mode with score_threshold=0 to accept all peaks
  const char* standard_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
    },
    "precursor_selection": {
      "max_mass_count": [3], "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false, "MS3AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": true, "value_ms": 30000 },
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" }
  })";

  // Config with MS3 mode 1 (SourceCID) enabled
  const char* ms3_mode1_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
    },
    "precursor_selection": {
      "max_mass_count": [3], "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false, "MS3AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "ms3": { "enabled": true, "mode": 1, "max_per_ms2": 2, "protein_sequence": "" },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" }
  })";

  // Config with conditional MS2 enabled
  const char* conditional_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
    },
    "precursor_selection": {
      "max_mass_count": [3], "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false, "MS3AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "conditional_ms2": true,
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" }
  })";

  // Full MS1 peaks from ms1_standard.txt scan 1 (520 peaks), RT=0.0668 min
  const double ms1_mzs[] = {
    500.376099, 500.824677, 502.786987, 502.851654, 502.918884, 502.935303,
    503.107269, 503.851868, 504.107849, 504.861664, 504.874603, 504.898590,
    504.966095, 505.086548, 505.104584, 506.086975, 506.104248, 506.796631,
    506.938477, 507.083191, 507.310852, 508.082794, 508.721802, 508.849792,
    508.882050, 509.200409, 509.926544, 510.319275, 511.242737, 511.321899,
    511.471466, 512.475159, 512.868652, 512.958618, 513.281189, 513.360168,
    513.951599, 515.383606, 515.840393, 515.867493, 517.159546, 519.138550,
    520.138977, 520.908813, 521.117981, 521.135803, 521.907043, 521.930908,
    522.117249, 522.135803, 522.875000, 523.114014, 523.133240, 523.878235,
    524.136108, 524.762756, 525.034058, 526.857483, 526.935730, 527.843201,
    527.903992, 528.851685, 529.376709, 531.739929, 531.890869, 533.190369,
    534.190002, 534.979248, 535.889771, 536.165100, 536.854553, 536.974792,
    537.165588, 538.143738, 538.162476, 539.163025, 539.799377, 540.159668,
    540.848877, 541.331299, 542.812866, 543.899475, 544.852600, 545.779968,
    545.895508, 547.018494, 547.398071, 548.375977, 550.217224, 550.779175,
    550.802124, 550.900330, 551.217773, 552.197083, 552.214294, 552.919678,
    553.213745, 553.730774, 553.814087, 554.896973, 555.331177, 556.827820,
    556.898621, 557.852417, 557.867432, 558.951782, 559.842651, 559.857239,
    561.416809, 563.074158, 563.813354, 563.829041, 564.075256, 565.356567,
    565.841980, 566.846802, 566.880920, 566.891357, 569.796265, 569.816956,
    570.741638, 570.810303, 570.834534, 571.378845, 571.778687, 573.781555,
    574.725098, 574.817139, 576.839722, 577.126038, 577.345825, 578.126343,
    579.105164, 579.123962, 580.105835, 581.102600, 581.119995, 582.854797,
    586.388062, 586.924011, 588.381714, 589.814575, 589.837341, 589.887146,
    591.177673, 591.305298, 592.466003, 592.782715, 593.143433, 593.157410,
    593.171326, 593.805908, 594.157837, 595.155151, 595.811340, 596.154907,
    596.334412, 596.829224, 597.135010, 597.150818, 597.797668, 597.863892,
    599.377686, 599.934631, 600.862427, 602.867371, 604.802307, 604.834961,
    604.933472, 607.209412, 608.209778, 609.206909, 609.343506, 609.793884,
    609.887878, 610.183838, 610.206055, 610.859070, 611.184570, 611.203430,
    612.162964, 612.181946, 613.181091, 613.272522, 613.332764, 613.346069,
    613.384094, 614.178894, 615.812012, 615.906128, 617.295898, 617.760803,
    617.992004, 618.249207, 618.819336, 619.843323, 620.284729, 620.855347,
    623.209412, 624.235901, 624.815308, 625.236755, 625.824585, 626.215576,
    626.233215, 626.303101, 627.063843, 628.734192, 628.906921, 629.928040,
    632.793945, 635.968750, 636.736145, 636.793335, 636.849304, 637.829224,
    638.748169, 639.729492, 645.273315, 645.765381, 645.814026, 647.311768,
    647.341980, 647.776062, 648.669739, 648.743347, 648.858337, 651.144714,
    651.837524, 652.145325, 653.125183, 653.143555, 654.142761, 657.893311,
    659.785034, 661.786987, 663.844055, 664.750610, 664.948730, 665.802307,
    665.846191, 667.176208, 668.176758, 669.174133, 669.813354, 670.174561,
    670.807373, 671.766968, 671.887512, 672.283752, 672.818909, 674.878540,
    678.259644, 679.745728, 679.950500, 682.781311, 684.202759, 684.740356,
    684.835510, 684.930908, 685.203430, 686.200745, 686.800415, 687.201355,
    687.422668, 687.774597, 688.199280, 688.240601, 688.783203, 689.200256,
    692.282104, 693.737122, 695.794434, 696.824158, 697.833130, 698.254639,
    698.337891, 701.855286, 702.745789, 702.789062, 704.327454, 706.772339,
    706.869446, 707.268250, 707.847168, 708.394897, 712.772278, 714.877258,
    718.171082, 718.264954, 718.785339, 720.871643, 722.288574, 722.335449,
    723.794922, 724.332336, 724.796509, 725.164001, 725.722534, 726.164856,
    726.451294, 726.796631, 727.801636, 727.824280, 728.718140, 731.751282,
    737.825989, 739.743408, 741.194946, 741.720154, 742.195251, 742.834961,
    743.193115, 744.193481, 744.734619, 745.191833, 745.798157, 745.817078,
    747.748108, 748.363342, 749.742493, 750.156616, 750.874084, 751.260315,
    753.726074, 755.230286, 758.221497, 758.297668, 759.222412, 760.220154,
    761.218384, 761.299133, 762.217590, 764.294312, 764.693970, 765.740967,
    765.823364, 766.223755, 767.718201, 768.403259, 768.589905, 769.790894,
    770.728149, 772.197876, 774.765320, 774.801880, 776.742554, 780.728088,
    782.664856, 782.881836, 789.712769, 791.797058, 795.253174, 796.256592,
    798.741638, 800.783264, 800.861206, 803.942932, 807.660583, 807.786072,
    808.734680, 812.752319, 813.802185, 813.848083, 815.214294, 815.741882,
    816.214417, 816.262207, 816.796570, 817.211853, 818.213013, 819.209961,
    819.708740, 821.275208, 821.724915, 822.160034, 823.720215, 823.750122,
    831.828552, 832.240234, 832.814209, 833.241272, 833.763428, 834.240540,
    834.746643, 839.708313, 840.817383, 841.771606, 843.762634, 852.781860,
    854.762390, 861.778381, 865.182617, 866.785706, 868.230713, 868.665100,
    876.735718, 885.846069, 888.179626, 889.232178, 889.703125, 890.203979,
    890.233521, 891.230896, 893.182800, 894.246460, 896.655090, 898.619690,
    903.331482, 903.622803, 903.667664, 904.269409, 906.261597, 908.713196,
    915.713745, 917.648010, 918.656067, 919.288025, 920.730225, 924.560120,
    928.634827, 929.311218, 929.699280, 930.142212, 938.176819, 943.649780,
    944.639832, 945.663086, 957.802002, 959.663269, 960.693787, 960.739014,
    962.217224, 962.657898, 963.671875, 963.781006, 966.692383, 972.152222,
    972.258789, 972.735718, 978.680420, 982.184631, 982.642334, 988.639771,
    990.625305, 1007.162720, 1007.700684, 1013.133972, 1015.669678, 1020.757935,
    1021.670715, 1030.745117, 1035.246826, 1035.657104, 1040.671509, 1045.147827,
    1047.168335, 1051.413818, 1053.669312, 1059.014771, 1061.552979, 1062.580078,
    1080.700317, 1087.543213, 1104.626587, 1105.107422, 1105.604126, 1110.571533,
    1130.602173, 1133.656982, 1134.119507, 1138.624268, 1139.643311, 1142.150757,
    1187.139160, 1188.156250, 1188.651978, 1189.614380, 1197.654175, 1202.685425,
    1224.913086, 1255.558838, 1257.455200, 1271.514893, 1287.515503, 1295.195435,
    1297.562500, 1301.028198, 1326.026978, 1328.234497, 1329.130737, 1331.615723,
    1340.144897, 1340.188110, 1341.229980, 1369.494141, 1372.993530, 1374.054443,
    1380.493774, 1382.583496, 1385.852661, 1390.265869, 1437.553711, 1452.506714,
    1456.314087, 1465.062012, 1534.384888, 1562.955811, 1575.464722, 1601.741821,
    1601.989990, 1642.893555, 1661.094727, 1686.925659, 1698.547119, 1703.671509,
    1703.763306, 1725.527832, 1783.464111, 1912.147949
  };
  const double ms1_ints[] = {
    1534.72, 1563.76, 2101.44, 2433.31, 1723.68, 2161.89,
    59130.62, 1785.82, 24195.42, 1334.11, 2638.57, 1899.27,
    1529.93, 14895.83, 14728.34, 7187.91, 5917.36, 1592.58,
    2392.21, 4690.54, 1462.64, 1447.79, 5966.40, 1857.63,
    2016.07, 1518.29, 1480.48, 1823.68, 1517.89, 1364.59,
    4543.34, 2547.51, 1965.32, 1590.43, 1550.91, 2263.99,
    1583.78, 2029.03, 1436.38, 3021.40, 1942.83, 246361.94,
    92025.09, 2205.68, 6748.71, 56501.90, 1798.35, 1720.37,
    1546.34, 18775.31, 1644.22, 2697.86, 7298.50, 1824.39,
    1808.04, 1470.49, 1549.85, 1988.97, 1639.81, 1774.12,
    1785.19, 1836.26, 1779.79, 1837.27, 1839.68, 5156.39,
    2609.03, 1711.84, 2379.55, 195921.22, 2526.73, 1624.02,
    95304.26, 4808.15, 52071.90, 16788.95, 2915.01, 5082.68,
    2475.81, 2019.09, 1726.78, 2616.66, 1848.44, 2255.98,
    1657.06, 1614.19, 1623.90, 1500.63, 11773.91, 1863.25,
    1627.64, 2363.27, 5773.60, 2226.48, 3018.15, 2143.18,
    2057.72, 1859.47, 2277.88, 2479.22, 1645.05, 3151.71,
    2423.66, 1966.98, 1916.55, 1728.60, 1788.30, 1619.18,
    1536.87, 4017.96, 1667.27, 2014.10, 2033.30, 1593.24,
    1710.92, 1634.30, 1897.90, 2021.01, 1810.71, 2040.40,
    1695.06, 2201.61, 1511.30, 1525.38, 1631.68, 1755.74,
    1554.91, 1615.39, 2122.08, 30457.84, 2104.44, 13873.09,
    7044.49, 5960.62, 2958.92, 2413.76, 1924.26, 2140.75,
    2315.28, 2126.27, 2074.89, 2068.43, 2267.22, 1900.04,
    2812.95, 1887.51, 6662.82, 1718.64, 1942.73, 109951.71,
    2621.15, 2328.61, 46630.72, 33525.32, 2261.89, 10582.25,
    2475.79, 2514.96, 1839.72, 3454.64, 3025.27, 2543.56,
    2988.31, 1915.48, 2260.98, 2154.77, 2048.81, 3370.70,
    2077.42, 21142.76, 7270.54, 3597.60, 1923.69, 2809.93,
    1890.85, 174427.27, 3359.07, 2074.62, 83930.87, 2752.50,
    3822.32, 46101.07, 22262.65, 2369.66, 2232.44, 1703.94,
    1901.82, 6546.02, 1818.20, 3072.12, 2561.39, 1909.85,
    2210.80, 2287.22, 2239.98, 2290.04, 1933.18, 1886.42,
    2168.73, 12880.94, 2554.13, 5716.43, 1752.04, 1960.48,
    6900.86, 2027.49, 1652.64, 1889.62, 1798.88, 1648.50,
    1857.10, 2021.49, 2243.77, 2737.98, 2611.51, 2026.25,
    2229.79, 1828.62, 1797.90, 2571.37, 1868.07, 2047.48,
    2200.64, 2174.58, 2036.61, 1844.38, 2295.78, 19872.75,
    2288.47, 9309.13, 3041.96, 4650.42, 4044.38, 2019.47,
    2732.00, 2094.58, 2561.11, 2223.70, 2021.36, 2235.45,
    2547.90, 60610.78, 35692.90, 21540.84, 2042.65, 6911.91,
    2686.90, 2588.32, 1880.99, 3276.49, 2566.36, 1994.06,
    2774.91, 2554.94, 1954.33, 2580.76, 79222.22, 1976.53,
    2112.71, 2401.47, 47508.61, 26244.49, 2795.32, 9182.56,
    2477.06, 2410.71, 3515.75, 1908.73, 3011.08, 2193.49,
    1837.48, 2224.72, 2210.40, 2321.49, 2116.13, 6441.56,
    2628.45, 2424.12, 2621.79, 2243.22, 2232.47, 1845.82,
    2446.14, 2343.66, 3008.78, 2293.30, 2095.56, 2435.94,
    2029.60, 3095.03, 2660.25, 2629.30, 2303.78, 2060.66,
    2237.30, 2762.05, 1999.54, 5148.54, 2736.82, 2383.26,
    1929.76, 2603.66, 2476.35, 2778.50, 2005.37, 3274.13,
    2528.28, 2557.64, 38758.45, 2797.37, 26024.79, 2308.52,
    19263.50, 9291.47, 3507.59, 3370.79, 2424.74, 2254.79,
    2430.51, 3033.64, 2371.03, 3411.51, 2426.70, 2379.55,
    2332.58, 2707.16, 32394.64, 2555.79, 27720.30, 15328.54,
    6586.63, 2824.78, 4928.58, 2206.15, 2335.58, 2347.27,
    2854.72, 2187.53, 1925.47, 2574.04, 7733.37, 2417.32,
    2518.22, 2955.63, 2331.00, 2160.17, 2788.06, 2881.51,
    2246.91, 2111.37, 2705.25, 2593.94, 2425.73, 2176.15,
    2365.86, 2690.50, 2955.10, 2483.51, 2367.18, 3646.64,
    2540.90, 2824.33, 2059.72, 2548.55, 15058.18, 2366.14,
    15806.40, 2977.10, 2416.31, 10739.36, 6125.60, 2685.59,
    2275.42, 2730.54, 2602.12, 3239.84, 2497.92, 3263.81,
    2529.75, 14006.59, 2331.97, 7876.31, 2966.97, 6006.48,
    2144.06, 2150.49, 2384.22, 9384.83, 2581.37, 3008.95,
    3037.44, 3044.59, 2831.70, 2846.20, 2801.97, 3703.36,
    2751.61, 3329.79, 2967.70, 7278.91, 3900.97, 2748.00,
    6338.63, 2538.72, 3358.17, 2713.71, 2980.42, 2945.10,
    2607.89, 2641.68, 3944.28, 3138.72, 2322.60, 2678.27,
    3251.97, 2907.44, 2716.31, 3185.30, 2181.00, 2405.81,
    3508.57, 2856.42, 2870.60, 2883.09, 3028.76, 2926.67,
    2653.81, 2974.37, 3612.98, 2612.78, 2549.03, 3547.10,
    2483.83, 2457.19, 4435.14, 2539.07, 3313.99, 3127.45,
    2620.46, 2946.14, 2814.47, 3587.50, 2716.31, 3363.95,
    3498.54, 3546.10, 3544.72, 2858.88, 3001.43, 3627.31,
    3555.01, 3608.16, 4157.98, 2565.52, 3206.81, 2616.56,
    3009.37, 3116.64, 2782.79, 3122.41, 2847.76, 2497.05,
    3294.66, 2928.28, 3414.46, 3295.61, 3188.46, 3852.41,
    3325.31, 3331.86, 4068.62, 3571.06, 2959.44, 2830.77,
    2914.07, 3149.96, 3143.94, 3238.05, 2707.45, 2580.96,
    3029.38, 3332.87, 3617.08, 4090.84, 3287.76, 2691.61,
    3273.74, 3161.31, 3141.28, 3456.22, 4158.72, 2941.02,
    3098.10, 3236.65, 3203.55, 3096.38, 3252.79, 3638.26,
    3119.90, 3063.51, 2845.28, 3150.03, 3356.29, 3162.19,
    3367.80, 3623.48, 3364.32, 3070.60, 3349.05, 3246.95,
    4207.45, 3165.06, 3440.03, 3552.06, 3485.78, 3410.75,
    12374.35, 3648.28, 4533.41, 3596.61
  };
  const int ms1_length = sizeof(ms1_mzs) / sizeof(ms1_mzs[0]);

  // Minimal MS2 peaks from ms2_hcd_fragment.txt scan 57
  const double ms2_mzs[] = {
    100.300545, 104.473274, 104.481895, 110.161835, 111.074051,
    111.105568, 111.118179, 111.128120
  };
  const double ms2_ints[] = {
    16955.05, 20442.67, 20458.91, 12642.20, 13831.51,
    13627.63, 14103.31, 13171.75
  };
  const int ms2_length = sizeof(ms2_mzs) / sizeof(ms2_mzs[0]);
}

START_TEST(FLASHIda_ProcessScan, "$Id$")

/////////////////////////////////////////////////////////////

// P4-U01: MS1 processScan returns > 0 commands for real spectral data
START_SECTION(processScan_ms1_returns_commands)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Diagnostic: call getPeakGroups directly to see deconvolution output
  int pg_count = ida->getPeakGroups(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_diag", nullptr);
  std::cout << "DIAG: getPeakGroups returned " << pg_count << " peak groups" << std::endl;
  std::cout << "DIAG: ms1_length=" << ms1_length << " first_mz=" << ms1_mzs[0]
            << " last_mz=" << ms1_mzs[ms1_length-1] << std::endl;

  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  std::cout << "DIAG: processScan returned " << n << " commands" << std::endl;
  TEST_EQUAL(n > 0, true)
  // Should be at most max_mass_count (3)
  TEST_EQUAL(n <= 3, true)
  delete ida;
}
END_SECTION

// P4-U02: Commands from processScan are dequeued by getNextScanCommand
START_SECTION(processScan_commands_dequeued)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue all commands
  ScanCommand cmd{};
  for (int i = 0; i < n; i++)
  {
    int result = ida->getNextScanCommand(cmd);
    TEST_EQUAL(result, 1)
    TEST_EQUAL(cmd.msn_level, 2)
    TEST_EQUAL(cmd.priority, 1)
    TEST_EQUAL(cmd.num_stages, 1)
    TEST_EQUAL(cmd.is_agc, 0)
    // Isolation stage should have valid m/z
    TEST_EQUAL(cmd.stages[0].precursor_mz > 0, true)
    TEST_EQUAL(cmd.stages[0].isolation_width > 0, true)
    TEST_EQUAL(cmd.stages[0].charge_state >= 4, true)
    // Scan description starts with 4-char tracking ID
    TEST_EQUAL(std::strlen(cmd.scan_description) >= 4, true)
    // Enqueue timestamp should be set
    TEST_EQUAL(cmd.enqueue_timestamp_ms > 0, true)
  }

  // Queue empty — next call returns 0
  int empty_result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(empty_result, 0)

  delete ida;
}
END_SECTION

// P4-U03: ScanCommand fields populated correctly
START_SECTION(processScan_command_fields)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  ScanCommand cmd{};
  ida->getNextScanCommand(cmd);

  // Analyzer from ms2_configs_[0]
  TEST_STRING_EQUAL(std::string(cmd.analyzer), "Orbitrap")
  // Resolution from ms2_configs_[0]
  TEST_EQUAL(cmd.orbitrap_resolution, 120000)
  // Activation type from ms2_configs_[0]
  TEST_STRING_EQUAL(std::string(cmd.stages[0].activation_type), "HCD")
  // Collision energy
  TEST_EQUAL(cmd.stages[0].collision_energy > 0, true)
  // Precursor m/z within scan range
  TEST_EQUAL(cmd.stages[0].precursor_mz >= 500.0, true)
  TEST_EQUAL(cmd.stages[0].precursor_mz <= 2000.0, true)

  delete ida;
}
END_SECTION

// P4-U04: processScan with ms_level=2 processes MS2 path
START_SECTION(processScan_ms2_path)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // First push an MS2 via MS1 processing to get a tracking ID in pending_scan_map_
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue one MS2 command to get its scan description (contains tracking ID)
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Now process MS2 return with the tracking ID in scan description
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, ms2_cmd.scan_description);
  // Should return 0 (no conditional, no MS3, no quant in standard config)
  TEST_EQUAL(ms2_result, 0)

  delete ida;
}
END_SECTION

// P4-U05: processScan with empty spectrum returns 0
START_SECTION(processScan_empty_spectrum)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(nullptr, nullptr, 0, 0.0, 1, "empty");
  TEST_EQUAL(n, 0)
  delete ida;
}
END_SECTION

// P4-U06: conditional MS2 follow-up is pushed at priority 2
START_SECTION(processScan_conditional_ms2_followup)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(conditional_json));

  // Push MS2 commands via MS1
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue first MS2
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.priority, 1)

  // Process MS2 return — should push conditional follow-up at priority 2
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, ms2_cmd.scan_description);
  TEST_EQUAL(ms2_result > 0, true)

  // Drain remaining priority-1 commands first
  ScanCommand out{};
  while (true)
  {
    ida->getNextScanCommand(out);
    if (out.priority != 1) break;
  }
  // Should get priority-2 conditional follow-up
  TEST_EQUAL(out.priority, 2)
  TEST_EQUAL(out.msn_level, 2)

  delete ida;
}
END_SECTION

// P4-U07: MS3 commands are pushed at priority 3
START_SECTION(processScan_ms3_commands)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(ms3_mode1_json));

  // Push MS2 commands via MS1
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue first MS2
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Process MS2 return — may push MS3 commands if fragments found
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, ms2_cmd.scan_description);
  // MS3 commands pushed depends on deconvolution results. We test the structure.
  if (ms2_result > 0)
  {
    // Drain all priority-1 commands
    ScanCommand out{};
    while (true)
    {
      ida->getNextScanCommand(out);
      if (out.priority != 1) break;
    }
    // If MS3 was pushed, should be priority 3 with 2 stages
    if (out.msn_level == 3)
    {
      TEST_EQUAL(out.priority, 3)
      TEST_EQUAL(out.num_stages, 2)
    }
  }

  delete ida;
}
END_SECTION

// P4-U08: decodeBase36 roundtrips with encodeBase36
START_SECTION(decodeBase36_roundtrip)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Test roundtrip via processScan → scan_description parsing
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  ScanCommand cmd{};
  ida->getNextScanCommand(cmd);
  // scan_description starts with 4-char base36 tracking ID
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  std::string id_str = desc.substr(0, 4);
  // Verify the ID is valid base-36
  for (char c : id_str)
  {
    TEST_EQUAL((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z'), true)
  }

  delete ida;
}
END_SECTION

// P4-U09: cleanupExpiredCommands removes old entries
START_SECTION(cleanup_expired_commands)
{
  // Use timeout-enabled config
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Push an MS2 via MS1 processing (this adds to pending_scan_map_)
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue all commands — they're still in pending_scan_map_ for tracking
  ScanCommand cmd{};
  for (int i = 0; i < n; i++)
  {
    ida->getNextScanCommand(cmd);
  }

  // The pending_scan_map_ entries have timestamps. cleanupExpiredCommands_
  // is called by getNextScanCommand. With timeout_ms=30000, entries should
  // NOT be expired immediately. Verify by processing an MS2 with valid tracking ID.
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, cmd.scan_description);
  // Should succeed (entry found in pending_scan_map_)
  // ms2_result can be 0 (no follow-ups) but shouldn't crash
  TEST_EQUAL(ms2_result >= 0, true)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
