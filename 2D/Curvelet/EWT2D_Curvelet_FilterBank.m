function mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,option)

%=======================================================================
% function mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,option)
% 
% This function generates the curvelet filter bank corresponding to the
% detected scales Bw and detected angles Bt.
%
% Inputs:
%   -Bw: list of detected scales (in ]0,pi[)
%   -Bt: list of detected angles (in [-3pi/4,pi/4])
%   -W: image width
%   -H: image height
%   -option: type of filterbank (see doc of curvelet transform)
%
% Outputs:
%   -mfb: cell containing the curvelet filter bank (mfb{1} is the 
%   lowpass filter. For option 1 and 2, the next mfb{s}{t} are the bandpass 
%   filters, s corresponds to the scales and t to the direction; while for 
%   option 3 we have mfb{t}{s})
%
% Author: Jerome Gilles
% Institution: SDSU - Department of Mathematics and Statistics
% Version: 1.0 - 2013
% Version: 2.0 - 2015
%=======================================================================

if rem(H,2)==0
    extH=1;
    H=H+1;
else
    extH=0;
end

if rem(W,2)==0
    extW=1;
    W=W+1;
else
    extW=0;
end

if option==1
    %% Independent detection of scales and angles
    
    % We compute gammaw
    gammaw=pi; 
    Npic=length(Bw);
    for k=1:Npic-1
        r=(Bw(k+1)-Bw(k))/(Bw(k+1)+Bw(k));
        if r<gammaw 
           gammaw=r;
        end
    end

    r=(pi-Bw(Npic))/(pi+Bw(Npic));
    if r<gammaw 
        gammaw=r; 
    end
    gammaw=(1-1/max(W,H))*gammaw;     %this ensure that gamma is chosen as strictly less than the min
    if gammaw>Bw(1)
        gammaw=(1-1/max(W,H))*Bw(1);
    end

    % We compute DTheta
    DTheta=2*pi; 
    Npic=length(Bt);
    for k=1:Npic-1
        r=(Bt(k+1)-Bt(k))/2;
        if r<DTheta 
           DTheta=r;
        end
    end

    r=(Bt(1)+pi-Bt(end))/2;
    if r<DTheta 
        DTheta=r; 
    end
    DTheta=(1-1/max(W,H))*DTheta; %this ensure that Dtheta is chosen as strictly less than the min

    % We prepare the grid of angles and radius
    [theta,radius]=EWT_CreateAngleGrid(W,H);

    % Build the filters
    mfb=cell(length(Bw)+1,1);
    
    % Generate the scaling function like for the Littlewood-Paley case
    mfb{1}=EWT2D_Curvelet_Scaling(radius,Bw(1),gammaw,W,H);

    % Generate each angular sector filter
    for s=1:length(Bw)-1
       mfb{s+1}=cell(length(Bt),1);
       for t=1:length(Bt)-1
            mfb{s+1}{t} = EWT_Angular_sector(theta,radius,Bt(t),Bt(t+1),Bw(s),Bw(s+1),gammaw,DTheta);
       end
       mfb{s+1}{length(Bt)} = EWT_Angular_sector(theta,radius,Bt(end),Bt(1)+pi,Bw(s),Bw(s+1),gammaw,DTheta);
    end

    % Generate the last ones
    mfb{length(Bw)+1}=cell(length(Bt),1);
    for t=1:length(Bt)-1
       mfb{length(Bw)+1}{t} = EWT_Angular_sector(theta,radius,Bt(t),Bt(t+1),Bw(end),2*pi,gammaw,DTheta);
    end
    mfb{length(Bw)+1}{length(Bt)} = EWT_Angular_sector(theta,radius,Bt(end),Bt(1)+pi,Bw(end),2*pi,gammaw,DTheta);
    
elseif option==2  
    %% Scales detected first and then angles per scales
    
    % We compute gammaw
    gammaw=1; 
    Npic=length(Bw);
    for k=1:Npic-1
        r=(Bw(k+1)-Bw(k))/(Bw(k+1)+Bw(k));
        if r<gammaw 
           gammaw=r;
        end
    end

    r=(pi-Bw(Npic))/(pi+Bw(Npic));
    if r<gammaw 
        gammaw=r; 
    end
    gammaw=(1-1/max(W,H))*gammaw; %this ensure that gamma is chosen as strictly less than the min


    % We compute DTheta
    DTheta=2*pi*ones(length(Bw),1);
    for s=1:length(Bw)
        Npic=length(Bt{s});
        for k=1:Npic-1
            r=(Bt{s}(k+1)-Bt{s}(k))/2;
            if r<DTheta(s)
               DTheta(s)=r;
            end
        end

        r=(Bt{s}(1)+pi-Bt{s}(end))/2;
        if r<DTheta(s) 
            DTheta(s)=r; 
        end
        DTheta(s)=(1-1/max(W,H))*DTheta(s); %this ensure that gamma is chosen as strictly less than the min
    end
 
    % We prepare the grid of angles and radius
    [theta,radius]=EWT_CreateAngleGrid(W,H);
    
    % Build the filters
    mfb=cell(length(Bw)+1,1);

    % Generate the scaling function
    mfb{1}=EWT2D_Curvelet_Scaling(radius,Bw(1),gammaw,W,H);

    % Generate each bandpass filter
    for s=1:length(Bw)-1
       mfb{s+1}=cell(length(Bt{s}),1);
       for t=1:length(Bt{s})-1
            mfb{s+1}{t} = EWT_Angular_sector(theta,radius,Bt{s}(t),Bt{s}(t+1),Bw(s),Bw(s+1),gammaw,DTheta(s));
       end
       mfb{s+1}{length(Bt{s})} = EWT_Angular_sector(theta,radius,Bt{s}(end),Bt{s}(1)+pi,Bw(s),Bw(s+1),gammaw,DTheta(s));
    end

    % Generate the last one
    mfb{length(Bw)+1}=cell(length(Bt{end}),1);
    for t=1:length(Bt{end})-1
       mfb{length(Bw)+1}{t} = EWT_Angular_sector(theta,radius,Bt{end}(t),Bt{end}(t+1),Bw(end),2*pi,gammaw,DTheta(end));
    end
    mfb{length(Bw)+1}{length(Bt{end})} = EWT_Angular_sector(theta,radius,Bt{end}(end),Bt{end}(1)+pi,Bw(end),2*pi,gammaw,DTheta(end));

else
    %% Angles detected first then scales per angles
    
    % We compute DTheta
    DTheta=2*pi; 
    Npic=length(Bt);
    for k=1:Npic-1
        r=(Bt(k+1)-Bt(k))/2;
        if r<DTheta 
           DTheta=r;
        end
    end

    r=(Bt(1)+pi-Bt(end))/2;
    if r<DTheta 
        DTheta=r; 
    end
    DTheta=(1-1/max(W,H))*DTheta; %this ensure that Dtheta is chosen as strictly less than the min
   
    % We compute gammaw
    gammaw=Bw{1}/2;
    for t=2:length(Bt)
        Npic=length(Bw{t});
        for k=1:Npic-1
            r=(Bw{t}(k+1)-Bw{t}(k))/(Bw{t}(k+1)+Bw{t}(k));
            if r<gammaw 
               gammaw=r;
            end
        end
        r=(pi-Bw{t}(end))/(pi+Bw{t}(end));
        if r<gammaw 
            gammaw=r; 
        end
    end
    gammaw=(1-1/max(W,H))*gammaw; %this ensure that gamma is chosen as strictly less than the min    
    
    % We prepare the grid of angles and radius
    [theta,radius]=EWT_CreateAngleGrid(W,H);

    % Build the filters
    mfb=cell(length(Bt)+1,1);

    % Generate the scaling function
    mfb{1}=EWT2D_Curvelet_Scaling(radius,Bw{1},gammaw,W,H);
    
    % Generate each angular sector filter
    for t=1:length(Bt)-1
       mfb{t+1}=cell(length(Bw{t+1})+1,1);
       % Generate the first one
       mfb{t+1}{1} = EWT_Angular_sector(theta,radius,Bt(t),Bt(t+1),Bw{1},Bw{t+1}(1),gammaw,DTheta);
       % Generate the other scales in the angular sector
       for s=1:length(Bw{t+1})-1
            mfb{t+1}{s+1} = EWT_Angular_sector(theta,radius,Bt(t),Bt(t+1),Bw{t+1}(s),Bw{t+1}(s+1),gammaw,DTheta);
       end
       mfb{t+1}{end} = EWT_Angular_sector(theta,radius,Bt(t),Bt(t+1),Bw{t+1}(end),2*pi,gammaw,DTheta);
    end

    % Generate the last one
    mfb{end}=cell(length(Bw{end})+1,1);
    % Generate the first one
    mfb{end}{1} = EWT_Angular_sector(theta,radius,Bt(end),Bt(1)+pi,Bw{1},Bw{end}(1),gammaw,DTheta);
    % Generate the other scales in the angular sector
    for s=1:length(Bw{end})-1
          mfb{end}{s+1} = EWT_Angular_sector(theta,radius,Bt(end),Bt(1)+pi,Bw{end}(s),Bw{end}(s+1),gammaw,DTheta);
    end
    mfb{end}{end} = EWT_Angular_sector(theta,radius,Bt(end),Bt(1)+pi,Bw{end}(end),2*pi,gammaw,DTheta);
end


% check if we need to cut back to the original size
if extH==1
    mfb{1}=mfb{1}(1:end-1,:);
    for s=2:length(mfb)
        for t=1:length(mfb{s})
            mfb{s}{t}=mfb{s}{t}(1:end-1,:);
        end
    end
end

if extW==1
    mfb{1}=mfb{1}(:,1:end-1);
    for s=2:length(mfb)
        for t=1:length(mfb{s})
            mfb{s}{t}=mfb{s}{t}(:,1:end-1);
        end
    end
end

% shift back the filters in the correct position for filtering
mfb{1}=ifftshift(mfb{1});
for s=2:length(mfb)
    for t=1:length(mfb{s})
        mfb{s}{t}=ifftshift(mfb{s}{t});
    end
end

% symmetrize the Nyquist frequencies for even size images
[H , W]=size(mfb{1});
if option<3
    if extH==1
       s=zeros(1,W);
       for t=1:length(mfb{end})
          mfb{end}{t}(floor(H/2)+1,2:floor(W/2))=mfb{end}{t}(floor(H/2)+1,2:floor(W/2))+mfb{end}{t}(floor(H/2)+1,end:-1:floor(W/2)+2);
          mfb{end}{t}(floor(H/2)+1,floor(W/2)+2:end)=mfb{end}{t}(floor(H/2)+1,floor(W/2):-1:2);
          s=s+mfb{end}{t}(floor(H/2)+1,:).^2;
       end
       %normalization to preserve the tight frame condition
       for t=1:length(mfb{end})
          mfb{end}{t}(floor(H/2)+1,2:floor(W/2))=mfb{end}{t}(floor(H/2)+1,2:floor(W/2))./sqrt(s(2:floor(W/2)));
          mfb{end}{t}(floor(H/2)+1,floor(W/2)+2:end)=mfb{end}{t}(floor(H/2)+1,floor(W/2)+2:end)./sqrt(s(floor(W/2)+2:end));
       end

    end

    if extW==1
       s=zeros(H,1);
       for t=1:length(mfb{end})
          mfb{end}{t}(2:floor(H/2),floor(W/2)+1)=mfb{end}{t}(2:floor(H/2),floor(W/2)+1)+mfb{end}{t}(end:-1:floor(H/2)+2,floor(W/2)+1);
          mfb{end}{t}(floor(H/2)+2:end,floor(W/2)+1)=mfb{end}{t}(floor(H/2):-1:2,floor(W/2)+1);
          s=s+mfb{end}{t}(:,floor(W/2)+1).^2;
       end
       %normalization to preserve the tight frame condition
       for t=1:length(mfb{end})
          mfb{end}{t}(2:floor(H/2),floor(W/2)+1)=mfb{end}{t}(2:floor(H/2),floor(W/2)+1)./sqrt(s(2:floor(H/2)));
          mfb{end}{t}(floor(H/2)+2:end,floor(W/2)+1)=mfb{end}{t}(floor(H/2)+2:end,floor(W/2)+1)./sqrt(s(floor(H/2)+2:end));
       end
    end
else
     if extH==1
       s=zeros(1,W);
       for t=2:length(mfb)
          mfb{t}{end}(floor(H/2)+1,2:floor(W/2))=mfb{t}{end}(floor(H/2)+1,2:floor(W/2))+mfb{t}{end}(floor(H/2)+1,end:-1:floor(W/2)+2);
          mfb{t}{end}(floor(H/2)+1,floor(W/2)+2:end)=mfb{t}{end}(floor(H/2)+1,floor(W/2):-1:2);
          s=s+mfb{t}{end}(floor(H/2)+1,:).^2;
       end
       %normalization to preserve the tight frame condition
       for t=2:length(mfb)
          mfb{t}{end}(floor(H/2)+1,2:floor(W/2))=mfb{t}{end}(floor(H/2)+1,2:floor(W/2))./sqrt(s(2:floor(W/2)));
          mfb{t}{end}(floor(H/2)+1,floor(W/2)+2:end)=mfb{t}{end}(floor(H/2)+1,floor(W/2)+2:end)./sqrt(s(floor(W/2)+2:end));
       end

    end

    if extW==1
       s=zeros(H,1);
       for t=2:length(mfb)
          mfb{t}{end}(2:floor(H/2),floor(W/2)+1)=mfb{t}{end}(2:floor(H/2),floor(W/2)+1)+mfb{t}{end}(end:-1:floor(H/2)+2,floor(W/2)+1);
          mfb{t}{end}(floor(H/2)+2:end,floor(W/2)+1)=mfb{t}{end}(floor(H/2):-1:2,floor(W/2)+1);
          s=s+mfb{t}{end}(:,floor(W/2)+1).^2;
       end
       %normalization to preserve the tight frame condition
       for t=2:length(mfb)
          mfb{t}{end}(2:floor(H/2),floor(W/2)+1)=mfb{t}{end}(2:floor(H/2),floor(W/2)+1)./sqrt(s(2:floor(H/2)));
          mfb{t}{end}(floor(H/2)+2:end,floor(W/2)+1)=mfb{t}{end}(floor(H/2)+2:end,floor(W/2)+1)./sqrt(s(floor(H/2)+2:end));
       end
    end
end