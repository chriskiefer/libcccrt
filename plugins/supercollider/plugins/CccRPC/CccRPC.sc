CccRPC : UGen {
	*ar { |input, highDim=10, lowDim=2|
		^this.multiNew('audio', input, highDim, lowDim);
	}
	checkInputs {
 		if (inputs.at(1).rate != 'scalar', {
 			^("highDim in not modulateable, please use a scalar value: " + inputs.at(1) + inputs.at(1).rate);
 		});
 		if (inputs.at(2).rate != 'scalar', {
 			^("lowDim in not modulateable, please use a scalar value: " + inputs.at(2) + inputs.at(2).rate);
 		});

		^this.checkValidInputs;
	}
}
