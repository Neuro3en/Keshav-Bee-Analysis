function y=lowpass(x,dt,Parameters)
% Return RC low-pass filter output samples (x), given input samples,
% time interval dt, and time constant RC
%
% Input:
%   x: the input of the low-pass filter, TxN matrix here T is number of
%   step along the time, N is the number of independent signals
%   dt: time step of the signal
%   Parameters: a structure
%   	-Parameters.RC: Time constant of the filter
    %Check input
    assert(ismatrix(x),'x need to be a 2D array')
    assert(size(x,1)>2,'x need to contain at least 2 rows')
    assert(isstruct(Parameters),'Parameters is not a structure')
    assert(isfield(Parameters,'RC'),'Parameters has no field name RC');
    %------
    y=NaN(size(x));
    Alpha=dt/(Parameters.RC+dt);
    %Set initial condition
    y(1,:)=x(1,:);
    %Low pass the signal with a 1st order low pass filter equation
    for t_i=2:size(x,1)
        y(t_i,:)=y(t_i-1,:)+Alpha*(x(t_i,:)-y(t_i-1,:));        
    end
end
