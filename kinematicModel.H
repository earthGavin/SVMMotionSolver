    // Motion is a vectorField of all moving boundary points
    vectorField motionBoundary(statPoints_.size(), vector::zero);

    scalar oldTime = mesh().time().value() - mesh().time().deltaT().value();
    scalar curTime = mesh().time().value();

    scalar alphaOld   = 0.0;
    scalar alphaCur   = 0.0;
    const scalar pi=3.141593;

    alphaOld = rotationAmplitude_*Foam::sin(2*pi*rotationFrequency_*oldTime);
    alphaCur = rotationAmplitude_*Foam::sin(2*pi*rotationFrequency_*curTime);

    vector translationVector
    (
        translationAmplitude_[0]*
        (
            Foam::sin(2*pi*translationFrequency_[0]*curTime)
          - Foam::sin(2*pi*translationFrequency_[0]*oldTime)
        ),
        translationAmplitude_[1]*
        (
            Foam::sin(2*pi*translationFrequency_[1]*curTime)
          - Foam::sin(2*pi*translationFrequency_[1]*oldTime)
        ),
        0
    );

    tensor RzOld
    (
        Foam::cos(alphaOld),  Foam::sin(alphaOld), 0,
       -Foam::sin(alphaOld),  Foam::cos(alphaOld), 0,
        0, 0, 1
    );

    tensor RzCur
    (
        Foam::cos(alphaCur),  Foam::sin(alphaCur), 0,
       -Foam::sin(alphaCur),  Foam::cos(alphaCur), 0,
        0, 0, 1
    );

    vectorField rotationField
    (
        (RzCur - RzOld)
      & (statPoints_ - initialRotationOrigin_)
    );

    motionBoundary = translationVector + rotationField;

    if (Pstream::parRun())
    {
	    localMotion_ = translationVector + ((RzCur - RzOld) & (localStatPoints_ - initialRotationOrigin_));
    }
