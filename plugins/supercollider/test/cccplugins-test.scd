{CccRPC.ar(SinOsc.ar(1000)).poll}.play


{CccRPC.ar(WhiteNoise.ar).poll}.play

{var w = XFade2.ar(WhiteNoise.ar, Saw.ar(100), MouseX.kr(-1,1)); CccRPC.ar(w, 10, 2).poll; w}.play