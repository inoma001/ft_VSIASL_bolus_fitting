function x = multiPID_vsiASL_bir8_timeseries_multi_global(ts,slice_no,M0)



% create a time vector for each tag-cntrl pair
xq = [550,550,550,550,660,660,660,660,800,800,800,800,960,960,960,960,1150,1150,1150,1150,1390,1390,1390,1390,1670,1670,1670,1670,2000,2000,2000,2000,2400,2400,2400,2400,2900,2900,2900,2900,3500,3500,3500,3500];
   
 
t=xq + 42.5*(slice_no-1) + 82; % addd on delay time (ms) for slice readout time and fixed delay from global crushing time to readout of the first slice
t_seconds=t/1000;

ts=ts'; %transpose timeseries vector

a_label=0.56;
a_BGS=0.95*0.95*0.95;
lambda=0.9;

t1b = 1.6; %T1 of blood in seconds
t1 = 1.3; %T1 of tissue in seconds


%objective function using Buxton multi compartment model with an extra macrovascular compartment - arterial + venous blood volume (aBV)

function fout = objectivefcn(x)
        
        bolus_duration = x(1); 

        f=x(2)*100; %cbf in ml/100g/min

        aBV=x(3);
        %limit aBV values to be positive
        if (aBV<0)
            aBV=0;
        end

        macro_duration = x(4)*1000; 

        bat=0;

        %upper limit of macro duration is bolus duration ( if aBV > 0)
        if (aBV>0) && ( (macro_duration/1000) > bolus_duration  )
            macro_duration =   bolus_duration*1000 ;
        end

        if (bolus_duration < xq(1)/1000 )
            bolus_duration = xq(1)/1000;
        end     
                
        tdash = xq ./1000;% label to crush time
        %correct for bat
        tdash=tdash-bat;
        tdash(le(tdash,0))=0;

        tdash(ge(tdash,bolus_duration)) = bolus_duration ; 

        % perfusion model as per Buxton 
        % f is flow per seconds 
        f_seconds = f / 6000;
        one_t1prime = 1/t1 + f_seconds/lambda;
        k = 1/t1b - one_t1prime;


        qp_t = ( exp(k.*t_seconds).*(1 - exp(-k.*tdash)) )./ (k.*tdash);

        asl_signal_mod = 2.*a_label.*a_BGS.*(M0/lambda).* f_seconds.*tdash.*(exp(-t_seconds./t1b)) .* qp_t;
        asl_signal_mod(isnan(asl_signal_mod))=0; %for when bat is greater than first LCT
        
       % create vector that is equal pure macro signal during macro bolus and zero when time is greater than the macro signal duration
        macro_contamination = t;
        macro_contamination(gt(macro_contamination,macro_duration)) = 0 ; 
        macro_vector_scaling =  (2 * a_label * a_BGS * M0 * exp(-t_seconds/t1b)) ./ 100 ; % / signal is boxcar that decays with T1 of blood
        macro_index = find(macro_contamination>0)  ;
        macro_contamination(gt(macro_contamination,0)) = macro_vector_scaling(macro_index) ;


        total_signal = ( aBV * macro_contamination ) + asl_signal_mod  ; %combination of signal models as per M Chappel et al. 2010 (perfusion is per tissue volume including aBV, CSF, GM, WM)

        fout =  sqrt (mean ( (ts - total_signal).^2 ) ) ; %rms error
       
end


%create vector of bolus macro bolus times for start points in minimisation
A = unique(t);
macro_bolus=A(1:end-1)./1000; 
%create multi-start vector for bolus duration, and macro duration as these start points have potentia to bias the fitting most
[X,Y,Z,V] = ndgrid(1.0:.75:3.25, 0.5, 0, macro_bolus);
 
W = [X(:),Y(:),Z(:),V(:)];
custpts = CustomStartPointSet(W);
x0 = [2.0, 0.5, 0.0, 0.0 ]; % bolus, cbf0, aBV, macro_duration

opts = optimoptions('fminunc','Display','off','Algorithm','quasi-newton'); % unconstrained minimisation of rms error (comparing ASL acquisitions with kinetic model)
problem = createOptimProblem('fminunc','objective', @objectivefcn,'x0',x0,'options',opts);

% Construct a MultiStart object 
ms = MultiStart;
ms.Display="final";
% minimise with multi-start
[x,~,~,~,solsms] = run(ms,problem,custpts);

% tidy up data

if (x(3) < 0) % don't let aBV be negative
    x(3) = 0;
end

%upper limit of macro duration is bolus duration ( if aBV > 0)
if (x(3)>0) && ( x(4) > x(1)  )
     x(4) =   x(1) ;
end

if (x(1) < xq(1)/1000) % if bolus is less first PID set equal to PID(1) mtahcing objective function
    x(1) = xq(1)/1000;
end

if (x(1) > xq(end)/1000) % if bolus greater than last PID set to equal last PID (as we are insensitive to bolus durations longer than last PID)
    x(1) = xq(end)/1000;
end

x(2)=x(2)*100; %cbf in ml/100g/min

% set bolus duration and cbf0 estimates to zero if cbf<0
if x(2)<=0
      x(1)=0;
      x(2)=0;
end

% set macro duration to zero if aBV=0 (as macro component only makes sense if aBV is greater than zero)
if x(3)==0
      x(4)=0;
end

%set minimum macro duration to 0
if x(4)<0
x(4)=0;
end


x

end
