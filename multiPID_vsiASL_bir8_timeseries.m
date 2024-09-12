function x = multiPID_vsiASL_bir8_timeseries(ts,slice_no,M0)



% create a time vector for each tag-cntrl pair
xq = [550,550,550,550,660,660,660,660,800,800,800,800,960,960,960,960,1150,1150,1150,1150,1390,1390,1390,1390,1670,1670,1670,1670,2000,2000,2000,2000,2400,2400,2400,2400,2900,2900,2900,2900,3500,3500,3500,3500];
   
t=xq + 42.5*(slice_no-1) + 82; % addd on delay time (ms) for slice readout time and fixed delay from global crushing time to readout of the first slice

ts=ts'; %transpose timeseries vector

t1_effective = 1600;
a_label=0.56;
a_BGS=0.95*0.95*0.95;
lambda=0.9;



%objective function using piecewise polynomial construction
    function f = objectivefcn(x)
        
        
        bolus_duration=x(1)*1000; % bolus duration in ms
        f=x(2)*100; %cbf in ml/100g/min
       
        
        breaks = [0 bolus_duration xq(end)]; 
        
        coefs = [1 0; 0 bolus_duration];
        pp = mkpp(breaks,coefs); 
   
        % create time axis that plateaus at the limit of the bolus
        % duration... this can then be used to generate the asl signal
        tdash = ppval(pp,xq)./1000; %bolus duration function in seconds - matches LCT time when LCT =< bolus duration, matches bolus duration when LCT > bolus duration
        
        asl_signal_mod = 2.*a_label.*a_BGS.*M0.*f.*tdash.*(exp(-t./t1_effective))./(6000*lambda);
        
        
        f =  sqrt (mean ( (ts - asl_signal_mod).^2 ) ) ; %rms error

    end



x0=[2.0, 0.5]; % initial guess for label duration is 2 seconds and cbf is 50ml/100g/min


options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton'); % unconstrained minimisation of rms error (comparing ASL acquisitions with kinetic model)
x = fminunc(@objectivefcn,x0,options);

x(2)=x(2)*100; %cbf in ml/100g/min


% remove estimates that are below lower bound
if x(1)<xq(1)/1000
      x(1)=0;
end

%limit max bolus duration to max LCT
if x(1)>xq(end)/1000
      x(1)=xq(end)/1000;
end
% set bolus duration and cbf0 estimates to zero if cbf<0
if x(2)<=0
      x(1)=0;
      x(2)=0;
end


end
