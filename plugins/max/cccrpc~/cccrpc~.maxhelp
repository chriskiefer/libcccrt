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
   620.0,
   520.0
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
      60.0
     ],
     "text": "Measures the dynamical complexity of a signal by projecting sliding windows through a fixed random matrix and counting the occupied cells of a low-dimensional histogram. Arguments: highDim (projection window, samples), lowDim (projection dimensions), maxWinSize (ms). Kiefer 2023, Sound and Music Computing."
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
     "text": "res 5",
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
      300.0,
      20.0
     ],
     "text": "attributes: analysis window (ms), hop (fraction of window),"
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
      285.0,
      300.0,
      20.0
     ],
     "text": "histogram resolution, projection hop (fraction of highDim)"
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
      320.0,
      110.0,
      22.0
     ],
     "text": "cccrpc~ 10 2 500",
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
      360.0,
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
      395.0,
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
      395.0,
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
      440.0,
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
      452.0,
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
      360.0,
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
      360.0,
      260.0,
      20.0
     ],
     "text": "float outlet: one value per analysis hop"
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
   }
  ],
  "dependency_cache": [],
  "autosave": 0
 }
}
