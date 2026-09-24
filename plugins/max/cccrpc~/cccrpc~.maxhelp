{
 "patcher": {
  "fileversion": 1,
  "appversion": {
   "major": 8,
   "minor": 5,
   "revision": 0,
   "architecture": "x64",
   "modernui": 1
  },
  "classnamespace": "box",
  "rect": [
   100.0,
   100.0,
   640.0,
   620.0
  ],
  "bglocked": 0,
  "openinpresentation": 0,
  "default_fontsize": 12.0,
  "default_fontface": 0,
  "default_fontname": "Arial",
  "gridonopen": 1,
  "gridsize": [
   15.0,
   15.0
  ],
  "gridsnaponopen": 1,
  "objectsnaponopen": 1,
  "statusbarvisible": 2,
  "toolbarvisible": 1,
  "boxes": [
   {
    "box": {
     "id": "obj-title",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      20.0,
      15.0,
      500.0,
      33.0
     ],
     "text": "cccrpc~  -  Random Projection Complexity",
     "fontsize": 20.0
    }
   },
   {
    "box": {
     "id": "obj-desc",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      20.0,
      50.0,
      560.0,
      75
     ],
     "text": "Measures the dynamical complexity of a signal by projecting sliding windows through a fixed random matrix and counting the occupied cells of a low-dimensional histogram. Kiefer 2023, Sound and Music Computing. Every setting is an attribute: type it in the box (@lowdim 4), send it as a message, or edit it in the inspector, where it is saved with the patch. The arguments of older patches (highDim lowDim maxWinSize maxLowDim) still work."
    }
   },
   {
    "box": {
     "id": "obj-noise",
     "maxclass": "newobj",
     "numinlets": 1,
     "numoutlets": 1,
     "patching_rect": [
      30.0,
      140.0,
      44.0,
      22.0
     ],
     "text": "noise~",
     "outlettype": [
      "signal"
     ]
    }
   },
   {
    "box": {
     "id": "obj-cycle",
     "maxclass": "newobj",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      100.0,
      140.0,
      66.0,
      22.0
     ],
     "text": "cycle~ 220",
     "outlettype": [
      "signal"
     ]
    }
   },
   {
    "box": {
     "id": "obj-c1",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      180.0,
      140.0,
      200.0,
      20.0
     ],
     "text": "noise = high complexity, sine = low"
    }
   },
   {
    "box": {
     "id": "obj-selnum",
     "maxclass": "number",
     "numinlets": 1,
     "numoutlets": 2,
     "patching_rect": [
      30.0,
      175.0,
      50.0,
      22.0
     ],
     "outlettype": [
      "",
      "bang"
     ],
     "minimum": 0,
     "maximum": 2
    }
   },
   {
    "box": {
     "id": "obj-c2",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      85.0,
      175.0,
      200.0,
      20.0
     ],
     "text": "1 = noise, 2 = sine, 0 = silence"
    }
   },
   {
    "box": {
     "id": "obj-sel",
     "maxclass": "newobj",
     "numinlets": 3,
     "numoutlets": 1,
     "patching_rect": [
      30.0,
      210.0,
      90.0,
      22.0
     ],
     "text": "selector~ 2 1",
     "outlettype": [
      "signal"
     ]
    }
   },
   {
    "box": {
     "id": "obj-m1",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      250.0,
      200.0,
      70.0,
      22.0
     ],
     "text": "winsize 25",
     "outlettype": [
      ""
     ]
    }
   },
   {
    "box": {
     "id": "obj-m2",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      330.0,
      200.0,
      70.0,
      22.0
     ],
     "text": "winsize 100",
     "outlettype": [
      ""
     ]
    }
   },
   {
    "box": {
     "id": "obj-m3",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      250.0,
      230.0,
      70.0,
      22.0
     ],
     "text": "hopsize 0.5",
     "outlettype": [
      ""
     ]
    }
   },
   {
    "box": {
     "id": "obj-m4",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      330.0,
      230.0,
      45.0,
      22.0
     ],
     "text": "res 10",
     "outlettype": [
      ""
     ]
    }
   },
   {
    "box": {
     "id": "obj-m5",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      385.0,
      230.0,
      50.0,
      22.0
     ],
     "text": "res 20",
     "outlettype": [
      ""
     ]
    }
   },
   {
    "box": {
     "id": "obj-m6",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      250.0,
      260.0,
      75.0,
      22.0
     ],
     "text": "rpchop 0.5",
     "outlettype": [
      ""
     ]
    }
   },
   {
    "box": {
     "id": "obj-c3",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      250.0,
      170.0,
      340.0,
      20.0
     ],
     "text": "attributes: send a message, or use the inspector"
    }
   },
   {
    "box": {
     "id": "obj-c4",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      250.0,
      378.0,
      370.0,
      48.0
     ],
     "text": "highdim, maxlowdim and maxwinsize set buffer sizes: changing them reallocates and restarts the analysis. lowdim above maxlowdim also reallocates."
    }
   },
   {
    "box": {
     "id": "obj-rpc",
     "maxclass": "newobj",
     "numinlets": 1,
     "numoutlets": 2,
     "patching_rect": [
      30.0,
      420.0,
      215.0,
      22.0
     ],
     "text": "cccrpc~ @highdim 10 @maxwinsize 500",
     "outlettype": [
      "signal",
      "float"
     ]
    }
   },
   {
    "box": {
     "id": "obj-snap",
     "maxclass": "newobj",
     "numinlets": 2,
     "numoutlets": 1,
     "patching_rect": [
      30.0,
      460.0,
      80.0,
      22.0
     ],
     "text": "snapshot~ 50",
     "outlettype": [
      "float"
     ]
    }
   },
   {
    "box": {
     "id": "obj-out",
     "maxclass": "flonum",
     "numinlets": 1,
     "numoutlets": 2,
     "patching_rect": [
      30.0,
      495.0,
      60.0,
      22.0
     ],
     "outlettype": [
      "",
      "bang"
     ]
    }
   },
   {
    "box": {
     "id": "obj-c5",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      100.0,
      495.0,
      180.0,
      20.0
     ],
     "text": "signal outlet via snapshot~"
    }
   },
   {
    "box": {
     "id": "obj-dac",
     "maxclass": "ezdac~",
     "numinlets": 2,
     "numoutlets": 0,
     "patching_rect": [
      30.0,
      540.0,
      45.0,
      45.0
     ]
    }
   },
   {
    "box": {
     "id": "obj-c6",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      85.0,
      552.0,
      200.0,
      20.0
     ],
     "text": "turn audio on"
    }
   },
   {
    "box": {
     "id": "obj-out2",
     "maxclass": "flonum",
     "numinlets": 1,
     "numoutlets": 2,
     "outlettype": [
      "",
      "bang"
     ],
     "patching_rect": [
      160.0,
      460.0,
      60.0,
      22.0
     ]
    }
   },
   {
    "box": {
     "id": "obj-c7",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      225.0,
      460.0,
      260.0,
      20.0
     ],
     "text": "float outlet: one value per analysis hop"
    }
   },
   {
    "box": {
     "id": "obj-m7",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      335.0,
      260.0,
      90.0,
      22.0
     ],
     "text": "downsample 1"
    }
   },
   {
    "box": {
     "id": "obj-m8",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      435.0,
      260.0,
      90.0,
      22.0
     ],
     "text": "downsample 8"
    }
   },
   {
    "box": {
     "id": "obj-m9",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      410.0,
      200.0,
      55.0,
      22.0
     ],
     "text": "lowdim 2"
    }
   },
   {
    "box": {
     "id": "obj-m10",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      475.0,
      200.0,
      55.0,
      22.0
     ],
     "text": "lowdim 4"
    }
   },
   {
    "box": {
     "id": "obj-n0",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      250.0,
      290.0,
      78.0,
      22.0
     ],
     "text": "normalize 0"
    }
   },
   {
    "box": {
     "id": "obj-n1",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      335.0,
      290.0,
      78.0,
      22.0
     ],
     "text": "normalize 1"
    }
   },
   {
    "box": {
     "id": "obj-n2",
     "maxclass": "message",
     "numinlets": 2,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "patching_rect": [
      420.0,
      290.0,
      78.0,
      22.0
     ],
     "text": "normalize 2"
    }
   },
   {
    "box": {
     "id": "obj-c8",
     "maxclass": "comment",
     "numinlets": 1,
     "numoutlets": 0,
     "patching_rect": [
      505.0,
      290.0,
      120.0,
      20.0
     ],
     "text": "raw / max / noise = 1"
    }
   },
   {
    "box": {
     "id": "obj-a1",
     "maxclass": "attrui",
     "attr": "highdim",
     "numinlets": 1,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "parameter_enable": 0,
     "patching_rect": [
      250.0,
      322.0,
      150.0,
      22.0
     ]
    }
   },
   {
    "box": {
     "id": "obj-a2",
     "maxclass": "attrui",
     "attr": "maxlowdim",
     "numinlets": 1,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "parameter_enable": 0,
     "patching_rect": [
      410.0,
      322.0,
      150.0,
      22.0
     ]
    }
   },
   {
    "box": {
     "id": "obj-a3",
     "maxclass": "attrui",
     "attr": "maxwinsize",
     "numinlets": 1,
     "numoutlets": 1,
     "outlettype": [
      ""
     ],
     "parameter_enable": 0,
     "patching_rect": [
      250.0,
      350.0,
      150.0,
      22.0
     ]
    }
   }
  ],
  "lines": [
   {
    "patchline": {
     "source": [
      "obj-noise",
      0
     ],
     "destination": [
      "obj-sel",
      1
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-cycle",
      0
     ],
     "destination": [
      "obj-sel",
      2
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-selnum",
      0
     ],
     "destination": [
      "obj-sel",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-sel",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m1",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m2",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m3",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m4",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m5",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m6",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-rpc",
      0
     ],
     "destination": [
      "obj-snap",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-snap",
      0
     ],
     "destination": [
      "obj-out",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-rpc",
      1
     ],
     "destination": [
      "obj-out2",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m7",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m8",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m9",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-m10",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-n0",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-n1",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-n2",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-a1",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-a2",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   },
   {
    "patchline": {
     "source": [
      "obj-a3",
      0
     ],
     "destination": [
      "obj-rpc",
      0
     ]
    }
   }
  ],
  "dependency_cache": [],
  "autosave": 0
 }
}
