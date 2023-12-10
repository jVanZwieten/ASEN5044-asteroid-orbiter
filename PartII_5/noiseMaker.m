function noise = noiseMaker(mu,R,nMeas)
    noise = mvnrnd(mu,R,nMeas);
end