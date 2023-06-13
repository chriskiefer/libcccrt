CccRPC : UGen {
	*ar { |input, highDim=10, lowDim=2, rpcHopSize=0.5, rpcRes=5, winSize=25, hopSize=0.5, maxWinSize=500|
		^this.multiNew('audio', input, highDim, lowDim, rpcHopSize, rpcRes, winSize, hopSize, maxWinSize);
	}
	checkInputs {
 		if (inputs.at(1).rate != 'scalar', {
 			^("highDim in not modulateable, please use a scalar value: " + inputs.at(1) + inputs.at(1).rate);
 		});
 		if (inputs.at(2).rate != 'scalar', {
 			^("lowDim in not modulateable, please use a scalar value: " + inputs.at(2) + inputs.at(2).rate);
 		});
 		if (inputs.at(7).rate != 'scalar', {
 			^("maxWinSize in not modulateable, please use a scalar value: " + inputs.at(7) + inputs.at(7).rate);
 		});

		^this.checkValidInputs;
	}
}
