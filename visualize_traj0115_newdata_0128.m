clear;
LW = 'linewidth';
prefix_arr = {...
    'him-5 male SP new data',...
};
tmpprefix = prefix_arr{1};
fname = ['analysis/',tmpprefix,'.mat'];
load(fname);

pos=[10 10 45 40];

mytransp = 0.2;
ms = 5;
myfs = 15;
mylw = 3;
fwidth = 1200; fheight = 2400;

for ii=1:Ndata
% for ii=3
    if(sum(ii == bad_inds))
        continue
    end
    ii
    
    x_stimulation = x{ii}(stimulation_frame_arr{ii});
    y_stimulation = y{ii}(stimulation_frame_arr{ii});

    trajid = trajids_arr{ii}
    figure('position',[400 400 3000 2200])%,'Visible','off');
    str_sgtitle=[tmpprefix, ' traj ',trajid];
    if(sum(ii == inds_not_find))
        str_sgtitle = [str_sgtitle,' not found'];
    end
    
    if(longer25_ids(ii)==1)
        disp('this one is longer than 25min')
        str_sgtitle = [str_sgtitle, ' (>25min)'];
    end
    
    if(longer15_ids(ii)==1)
        disp('this one is longer than 15min')
        str_sgtitle = [str_sgtitle, ' (>15min)'];
    end
    
    cur_set = setids(ii);
    target_pos = target_pos_set{cur_set};
    theta = linspace(0, 2*pi, 200);
    
    % set(gcf,'position',pos)
    sgtitle(str_sgtitle,'fontweight','normal','fontsize',myfs*2,'interpreter','none')
    
    % color by time
    subplot(5,4,1)
    % figure('Position',[200 300 600 500])
    scatter(x{ii},y{ii},ms,time{ii},'filled','MarkerFaceAlpha',mytransp);
    htime = gca;
    hold on
    scatter(x_stimulation, y_stimulation, ms*10, 'k','o',LW,mylw);
    scatter(target_pos(1),target_pos(2),100,'rs',LW,2)
    circleX = target_pos(1) + target_pos(3)*cos(theta);
    circleY = target_pos(2) + target_pos(3)*sin(theta);
    plot(circleX, circleY, 'r-');
    box on
    set(gca,'fontsize',myfs)
    cb = colorbar();
    clim([0, 22]);
    ylabel(cb, 'time [min]')
    axis('equal')
    
    % color by distance
    subplot(5,4,5)
    scatter(x{ii},y{ii},ms,distance{ii},'filled','MarkerFaceAlpha',mytransp)
    hdist = gca;
    hold on
    scatter(x_stimulation, y_stimulation, ms*10, 'k','o',LW,mylw);
    scatter(target_pos(1),target_pos(2),100,'rs',LW,2)
    circleX = target_pos(1) + target_pos(3)*cos(theta);
    circleY = target_pos(2) + target_pos(3)*sin(theta);
    plot(circleX, circleY, 'r-');
    box on
    set(gca,'fontsize',myfs)
    cb = colorbar();
    clim([0, 20]);
    ylabel(cb,'distance')
    axis('equal')

    subplot(5,4,[6,7,8])
    tt=time{ii};
    plot(tt, distance{ii}, LW,mylw)
    hold on 
    scatter(tt(stimulation_frame_arr{ii}), distance{ii}(stimulation_frame_arr{ii}), ms*10,'k','o',LW,mylw);
	title('distance','fontweight','normal')
    xlim([0,tt(end)])
    % xlabel('time [min]')
    ylabel('R')
    set(gca,'fontsize',myfs)
    box on


    % color by speed

    subplot(5,4,9)
    vv = v{ii};
    vv(vv>50)=0;
    % vmin = min(vv);
    % vmax = max(vv);
    % vdiff = vmax - vmin;
    % vmin = vmin - vdiff*0.1;
    % vmax = vmax + vdiff*0.1;
    vmin = 0;
    vmax = 15;
    scatter(x_mid{ii},y_mid{ii},ms,vv,'filled','MarkerFaceAlpha',mytransp)
    hspeed = gca;
    hold on
    scatter(x_stimulation, y_stimulation, ms*10, 'k','o',LW,mylw);
    scatter(target_pos(1),target_pos(2),100,'rs',LW,2)
    circleX = target_pos(1) + target_pos(3)*cos(theta);
    circleY = target_pos(2) + target_pos(3)*sin(theta);
    plot(circleX, circleY, 'r-');
    box on
    set(gca,'fontsize',myfs)
    cb = colorbar();
    clim([vmin, vmax]);
    ylabel(cb,'speed')
    axis('equal')

    subplot(5,4,[10,11,12])
    tt=time{ii};
    plot(time_mid{ii}, vv, LW,mylw)
    hold on 
    scatter(tt(stimulation_frame_arr{ii}), vv(stimulation_frame_arr{ii}), ms*10,'k','o',LW,mylw);
	title('speed','fontweight','normal')
    xlim([0,tt(end)])
    % xlabel('time [min]')
    ylabel('v')
    set(gca,'fontsize',myfs)
    box on

    % color by straightness

    subplot(5,4,13)
    scatter(x_mid{ii},y_mid{ii},ms,straightness{ii},'filled','MarkerFaceAlpha',mytransp)
    hold on
    scatter(x_stimulation, y_stimulation, ms*10, 'k','o',LW,mylw);
    hstraight = gca;
    scatter(target_pos(1),target_pos(2),100,'rs',LW,2)
    circleX = target_pos(1) + target_pos(3)*cos(theta);
    circleY = target_pos(2) + target_pos(3)*sin(theta);
    plot(circleX, circleY, 'r-');
    box on
    set(gca,'fontsize',myfs)
    cb = colorbar();
    clim([0 40]);
    ylabel(cb,'1/straightness')
    axis('equal')

    subplot(5,4,[14,15,16])
    tt=time{ii};
    plot(time_mid{ii}, straightness{ii}, LW,mylw)
    hold on 
    scatter(tt(stimulation_frame_arr{ii}), straightness{ii}(stimulation_frame_arr{ii}), ms*10,'k','o',LW,mylw);
	title('1/straightness','fontweight','normal')
    % xlabel('time [min]')
    xlim([0,tt(end)])
    ylabel('D_R')
    set(gca,'fontsize',myfs)
    box on


    % color by direction correctness


    subplot(5,4,17)
    scatter(x_mid{ii},y_mid{ii},ms,direction_correct{ii},'filled','MarkerFaceAlpha',mytransp)
    hold on
    scatter(x_stimulation, y_stimulation, ms*10, 'k','o',LW,mylw);
    hdirection = gca;
    colormap(gca,bluewhitered)
    clim([-1, 1])
    scatter(target_pos(1),target_pos(2),100,'rs',LW,2)
    circleX = target_pos(1) + target_pos(3)*cos(theta);
    circleY = target_pos(2) + target_pos(3)*sin(theta);
    plot(circleX, circleY, 'r-');
    box on
    set(gca,'fontsize',myfs)
    cb = colorbar();
    ylabel(cb, 'cos(\theta)')
    axis('equal')


    subplot(5,4,[18,19,20])
    tt=time_mid{ii};
    plot(tt, direction_correct{ii}, LW,mylw);
    hold on
    plot(tt, zeros(1,length(time_mid{ii})), 'k', LW,mylw);
    scatter(tt(stimulation_frame_arr{ii}), 0, ms*10,'k','o',LW,mylw);
	title('direction','fontweight','normal')
    xlabel('time [min]')
    xlim([0,tt(end)])
    ylabel('cos\theta')
    set(gca,'fontsize',myfs)
    box on

    colormap(htime, parula);
    colormap(hdist, hsv);
    colormap(hspeed, turbo);
    colormap(hstraight, flipud(winter));
%     colormap(hdirection, bluewhitered);

    figname=[prefix,'/',prefix,trajid,'.png'];
%     saveas(gcf, figname);
    print(gcf, figname, '-dpng', '-r300');
    
%     pause;
    close all;
end