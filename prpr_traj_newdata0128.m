%% load the data excel, already texttocolumned by macros.
%  then process it.

clear;

path = 'analysis/'
prefix = 'him-5 male SP new data';
tmpprefix = prefix;

newrun = 0;
if newrun
    mkdir([prefix]);
    fname = [prefix,'.xlsx'];
    data=readmatrix(fname,'numheaderlines',7);
    header_ids=readcell(fname,'range','A6:GX6');
    find_target = readcell(fname, 'range','A1:GX1');
    stimulation_frame = readcell(fname, 'Range','A2:GX2');
    find_target_frame = readcell(fname, 'Range','A3:GX3');
    data(isnan(data))=0;
    Ndata = ((size(data,2)-1)/2);
else
    load([path,prefix,'.mat'])
end

% this offset is for computing Dq
offset = 40;
offset2 = offset*2;
% % % % % % % % % % % % % % % % % 
path = 'analysis/';
Nframe=40;
polyorder = 3;
framesize = 95;
bad_inds = [2, 53, 54, 79, 20, 100, 65, 13,...
    44, 51, 68, 63, 34, 76, 85, 43]; % just remove these data


good_inds = [];
nkeep = 1;

dt_sec = 0.1333333333;

x = cell(1,Ndata); y = cell(1,Ndata);
x_mid = cell(1,Ndata); y_mid = cell(1,Ndata);
x_Dq=cell(1,Ndata); y_Dq=cell(1,Ndata);
lendata = size(data,1);
t = dt_sec*(1:lendata)/60; 
NN=zeros(1,Ndata);
v=cell(1,Ndata);
time=cell(1,Ndata);
time_mid=cell(1,Ndata);
time_Dq=cell(1,Ndata);
msdback=cell(1,Ndata);
drdts_arr=cell(1,Ndata);
v_raw=cell(1,Ndata);
drdts_arr_raw=cell(1,Ndata);

vel = zeros(lendata-Nframe+1,Ndata);
trajids_arr = cell(1,Ndata);
stimulation_frame_arr = cell(1,Ndata);
find_target_frame_arr = cell(1,Ndata);
stimulation_time_arr = cell(1,Ndata);
find_target_time_arr = cell(1,Ndata);
search_duration = cell(1,Ndata);
distance = cell(1,Ndata);
speed = cell(1,Ndata);
straightness = cell(1,Ndata);
direction_correct = cell(1,Ndata);

nkeep = 1;
longer15_ids = zeros(1,Ndata);
longer25_ids = zeros(1,Ndata);
dist_stim_target_arr = zeros(1,Ndata);

Nsets = 23;
inds_not_find = [];
jj = 1;

setids = zeros(1,Ndata);
setinds = cell(1,Nsets);
for jj=1:Nsets
    setinds{jj} = [];
end
% first get target position of each set
for ii=1:Ndata
    trajid = header_ids{2*ii};
    pattern = '\d+-';
    matchset = regexp(trajid, pattern, 'match');
    matchset = matchset{1};
    cur_set = str2num(matchset(1:end-1));
    setids(ii) = cur_set;
    setinds{cur_set} = [setinds{cur_set}; ii];
end

target_pos_set = cell(1, Nsets);
for jj=1:Nsets
    this_set = setinds{jj};
    ndata = length(this_set);
    if(ndata==0)
        continue
    end
    targetx_arr = [];
    targety_arr = [];
    for kk=1:ndata
        ii = this_set(kk)
        indx = 2*ii;
        indy = 2*ii+1;
        xx = data(:,indx)./1000;
        yy = data(:,indy)./1000;
        if(is_empty_data(xx))
            continue
        end
        [tt,xx,yy]=trimtrailingzeros(t,xx,yy);
        [xx,yy] = fillgap(tt,xx,yy);
        %%%%%% now process the header infos
        stim_frame = stimulation_frame{2*ii};
        target_frame = find_target_frame{2*ii};
        if(target_frame=='N')
            continue
        end
        if(target_frame>length(xx))
            target_frame = length(xx);
        end
        targetx_arr = [targetx_arr; xx(target_frame)];
        targety_arr = [targety_arr; yy(target_frame)];
    end
    scatter(targetx_arr, targety_arr, 100,'g','filled')
    hold on
    npoints = length(targetx_arr);
    if(npoints==0)
        bad_inds = [bad_inds, this_set(:)'];
        continue
    end
    if(npoints<3)
        disp(['only ',num2str(npoints),' points!']);
        xc = mean(targetx_arr);
        yc = mean(targety_arr);
        r = sqrt((xc-targetx_arr(1))^2+(yc-targety_arr(1))^2);
    end
    if(npoints==3)
        [xc, yc, r] =  get_center_3points(targetx_arr, targety_arr);
    end
    if(npoints>3)
        [xc, yc, r] = get_center(targetx_arr, targety_arr);
    end
    if(r>6)
        xc = median(targetx_arr);
        yc = median(targety_arr);
        dist = (targetx_arr-xc).^2+(targety_arr-yc).^2;
        sqrt(dist)
        r = mean(dist);
        if(r>6)
            r = median(sqrt(dist));
        end
    end
    disp(['r=',num2str(r)]);
    % this_set
    hold on
    scatter(xc,yc,100,'ro','filled')
    theta = linspace(0, 2*pi, 200);
    circleX = xc + r*cos(theta);
    circleY = yc + r*sin(theta);
    plot(circleX, circleY, 'b-');
    axis equal
    % pause
    close all
    target_pos_set{jj} = [xc,yc,r];
end

% return
% bad_inds

% now we have the good target positions
% 
inds_not_find = [];
jj=1;
for ii=1:Ndata
% for ii=37
    % disp(ii)
    cur_set = setids(ii);
    trajids_arr{ii} = header_ids{2*ii};

    if(sum(ii==bad_inds))
        disp(ii)
        continue
    end
    indx = 2*ii;
    indy = 2*ii+1;
    xx = data(:,indx)./1000;
    yy = data(:,indy)./1000;
    if(is_empty_data(xx))
        bad_inds = [bad_inds,ii];
        continue
    end

    %%%%%% now process the header infos
    stim_frame = stimulation_frame{2*ii};
    target_frame = find_target_frame{2*ii};
    target_pos = target_pos_set{cur_set};

    [tt,xx,yy]=trimtrailingzeros(t,xx,yy);
    [xx,yy] = fillgap(tt,xx,yy);
    
    NN(ii)=length(tt);
    if(target_frame=='N')
        % continue
        inds_not_find(jj)=ii;
        jj = jj+1;
        target_frame = NN(ii);
    end
    if(target_frame>NN(ii))
        disp(['data ',num2str(ii),' has target_frame=',num2str(target_frame),' but only lenseries=',num2str(NN(ii))]);
        disp('set it to lenseries');
        target_frame = NN(ii);
    end

    find_target_frame_arr{ii} = target_frame;
    stimulation_frame_arr{ii} = stim_frame;
    stimulation_time_arr{ii} = t(stimulation_frame_arr{ii});
    find_target_time_arr{ii} = t(find_target_frame_arr{ii});
    search_duration{ii} = t(find_target_frame_arr{ii}) - t(stimulation_frame_arr{ii});
    if(search_duration{ii}>25)
        longer25_ids(ii) = 1;
    elseif(search_duration{ii}>15)
        longer15_ids(ii) = 1;
    end

    if(stim_frame > NN(ii) || stim_frame>target_frame)
        bad_inds = [bad_inds, ii];
        disp([num2str(ii),' has stim_frame=',num2str(stim_frame),' but lenseries=', num2str(NN(ii)), ' and target_frame=',num2str(target_frame)]);
        disp('discard this data')
        continue
    end

    % if(ii==18)
    %     disp(ii)
    % end
    dist_stim_target = sqrt((xx(stim_frame)-target_pos(1))^2+...
        (yy(stim_frame)-target_pos(2))^2);
    dist_stim_target_arr(ii) = dist_stim_target;
    
    stim_frame_range = stim_frame-50:stim_frame+20;

    % now do Savitzky-Golay filter to smooth out the wiggles
    x_smooth = sgolayfilt(xx, polyorder, framesize);
    y_smooth = sgolayfilt(yy, polyorder, framesize);

    % size(xx)
    % size(x_smooth)
    % 
    % plot(xx,yy,LW,2)
    % hold on
    % % plot(x_smooth, y_smooth, LW,2)
    % scatter(x_smooth,y_smooth,ms,tt,'filled','MarkerFaceAlpha',mytransp);
    xx = x_smooth;
    yy = y_smooth;

    x_mid{ii} = (xx(1:end-1)+xx(2:end))./2;
    y_mid{ii} = (yy(1:end-1)+yy(2:end))./2;
    x{ii} = xx;
    y{ii} = yy;

    time_mid{ii}=(t(1:NN(ii)-1)+t(2:NN(ii)))./2;
    time{ii} = t(1:NN(ii));
    
    dt=t(2)-t(1);
    dx=xx(2:end)-xx(1:end-1);
    dy=yy(2:end)-yy(1:end-1);
    ds=sqrt(dx.^2+dy.^2);
    tmpv = ds(:)./dt(:);
    tmpv(tmpv<1e-4)=1e-4;
    tmpv(stim_frame_range) = (tmpv(stim_frame_range(1)-1)+tmpv(stim_frame_range(end)+1))/2;
    tmpvmovmean = movmean(tmpv,Nframe);
    v{ii}=tmpvmovmean;
    v_raw{ii}=tmpv;
    lenii = length(tmpv);
    
    Rx = target_pos(1) - x_mid{ii};
    Ry = target_pos(2) - y_mid{ii};
    R = sqrt(Rx.^2+Ry.^2);
    vx = dx ./ dt;
    vy = dy ./ dt;
    vmag = tmpv;
    costheta = (vx .* Rx + vy .* Ry) ./ vmag ./ R;
    costheta(end)=1;
    costhetamovmean = movmean(costheta,Nframe);
    direction_correct{ii} = costheta;
    
    vx = movmean(vx,offset2);
    vy = movmean(vy,offset2);
    direction_absolute = atan2(vy, vx);
    theta = angle(vx + vy*1i);
    dtheta = abs(theta(2+offset2:end) - theta(1:end-1-offset2));
    dtheta = min(dtheta, 2*pi-dtheta);
    % dtheta_dt = (dtheta)./(dt*(offset2+1));
    Dq = dtheta.^2./(dt*(offset2+1));
    
    % x_Dq{ii} = x{ii}(2+offset:end-1-offset);
    % y_Dq{ii} = y{ii}(2+offset:end-1-offset);
    time_Dq{ii} = time{ii}(2+offset:end-1-offset);
    Dq = interp1(time_Dq{ii},Dq, time_mid{ii})';
    
    straightness{ii} = Dq;
    direction_absolute_arr{ii} = direction_absolute;

    vel(1:lenii,ii) = v{ii}./max(v{ii});
    msd=(xx - target_pos(1)).^2 + (yy - target_pos(2)).^2;
    dd = sqrt(msd);
    distance{ii} = movmean(dd, Nframe);

    dr = (sqrt(msd(2:end))-sqrt(msd(1:end-1)));
    drdts=dr(:)./dt(:);
    drdts(stim_frame_range) = (drdts(stim_frame_range(1)-1)+drdts(stim_frame_range(end)+1))/2;
    drdts_arr{ii} = movmean(drdts, Nframe);
    drdts_arr_raw{ii} = drdts;
end
good_inds = setdiff(1:Ndata,bad_inds);

save([path,prefix,'.mat']);
% prefix = [prefix,'/',prefix];
% visualize_traj0115_newdata

function [xc, yc, r] = get_center_3points(x,y)
    midpoint1 = [(x(1) + x(2))/2, (y(1) + y(2))/2];
    slope1 = (y(2) - y(1)) / (x(2) - x(1));    
    midpoint2 = [(x(2) + x(3))/2, (y(2) + y(3))/2];
    slope2 = (y(3) - y(2)) / (x(3) - x(2));
    perpSlope1 = -1 / slope1;
    perpSlope2 = -1 / slope2;
    b1 = midpoint1(2) - perpSlope1 * midpoint1(1);
    b2 = midpoint2(2) - perpSlope2 * midpoint2(1);
    A = [-perpSlope1, 1; -perpSlope2, 1];
    B = [b1; b2];
    center = A \ B;
    xc = center(1); yc = center(2);
    r = sqrt((x(1) - xc)^2 + (y(1) - yc)^2);
end

function [xc, yc, r] = get_center(x, y)
    Npoints = length(x);
    N = nchoosek(Npoints, 3);
    xc = zeros(1,N);
    yc = zeros(1,N);
    r = zeros(1,N);
    xcombo = nchoosek(x, 3);
    ycombo = nchoosek(y, 3);
    for ii=1:N
        [xc(ii),yc(ii),r(ii)] = get_center_3points(xcombo(ii,:),ycombo(ii,:));
        scatter(xc(ii),yc(ii),100,'bo','filled')
        hold on
        theta = linspace(0, 2*pi, 200);
        circleX = xc(ii) + r(ii)*cos(theta);
        circleY = yc(ii) + r(ii)*sin(theta);
        plot(circleX, circleY, 'k-');
    end
    xc = median(xc);
    yc = median(yc);
    r = median(r);
end

function tt = is_empty_data(xx)
    a=(xx==0);
    tt = (sum(a)==length(a));
end
