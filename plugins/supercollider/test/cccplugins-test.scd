{CccRPC.ar(SinOsc.ar(1000)).poll}.play


{CccRPC.ar(WhiteNoise.ar).poll}.play

(
{
	var w = XFade2.ar(WhiteNoise.ar, LFPulse.ar(100), MouseX.kr(-1,1));
	var rpc=CccRPC.ar(w, 8, 4, 0.5, 3, 10,0.5).poll;
	[w, rpc/RunningMax.ar(rpc)]
}.scope
)


(
{
	CccRPC.ar(SoundIn.ar(0), 8, 4, 0.5, 3, 300,0.5).poll;
	0
}.play
)


(
Gui
)