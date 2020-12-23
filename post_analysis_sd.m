function post_analysis_sd

% ** Run the PanicPackage code first!
% run these functions in the order they are listed

% coarse_sd3(); % 1 coarse the dataset so that we sample at 0.01 s
% get_evac_times(); % 2 extract evacuation times for 90% of the crowd
% replay_a_trial(); % for verification 
% plot_traj(); % for verification 
% calculate_distmat(); % 3 this takes time! go for a coffee break!
% calculate_exposure_from_IA(16); % 4 takes the tau values [1 5 8 10 16];
composite_plot(); % 5 draw the composite plot for paper


function coarse_sd3()
% subsample the trajectories at 0.01 s intervals

files=dir('../PanicPackage/dump/sd_u*100agents_*.dat3');
dt1=1/10; % 1/100;



for ff=1:size(files,1)
    data=dlmread(['../PanicPackage/dump/', files(ff).name]);
    fprintf('processing %s ...\n', files(ff).name);
    % sort time into 1/100 dt increments

    simtime=unique(data(:,1));

    keeptime=simtime(1);
    for ii=2:numel(simtime)
        simtime1=simtime(simtime-keeptime(end)-dt1>0);
        [val, idx]=min(simtime1-keeptime(end)-dt1);
        keeptime=[keeptime; simtime1(idx)];
    end

    data2=[];
    for ii=1:numel(keeptime)
        data2=[data2; data(data(:,1)==keeptime(ii),:)];
    end
    
    % remove any 0,0 points
    data2=data2(data2(:,3)~=0,:);
    
    % write to csv
    csvwrite(['../PanicPackage/dump/', files(ff).name(1:end-4), 'csv'], data2);
end

function get_evac_times()

% number of remaining agents
nr=2;

files=dir('../PanicPackage/dump/sd_u*100agents_*.dat');

if contains(files(1).name, '100agents')
    nr=10;
end

for ff=1:size(files,1)
    data=dlmread(['../PanicPackage/dump/', files(ff).name],'\t', 1,0);
    fprintf('processing %s ...\n', files(ff).name);
    % sort time into 1/100 dt increments

    idx=find(data(:,3)<=nr);
    evactime=data(idx(1),2);
    % write to csv
    csvwrite(['../PanicPackage/dump/evactime_', files(ff).name(1:end-3), 'csv'], evactime);
end


function replay_a_trial

th=-pi:.1:pi;

files=dir('./sd_seed*.csv');
ff=1;



traj=csvread(['./', files(ff).name]);
ts=unique(traj(:,1));
ids=unique(traj(:,2));
siz=0.5+rand(1,numel(ids))*.2;
for k=1:numel(ts)
    figure(1); gcf; clf;
    pos=traj(traj(:,1)==ts(k),3:4);
    for jj=1:numel(ids)
        fill(pos(jj,1)+siz(jj).*cos(th)/2, pos(jj,2)+siz(jj).*sin(th)/2, ...
            0, 'LineStyle','none');
        hold on;
    end
    drawnow
end
axis([0 15 0 10]);
box on


function plot_traj

files=dir('../PanicPackage/dump/sd_u*.csv');
figure(1); gcf; clf;

for ff=1:size(files,1)
    traj=csvread(['../PanicPackage/dump/', files(ff).name]);
    ids=unique(traj(:,2));
    subplot(4,size(files,1)/4,ff);
    for ii=1:numel(ids)
        traj1=traj(traj(:,2)==ids(ii), :);
        [val, idx]=sort(traj1(:,1));
        plot(traj1(idx,3), traj1(idx,4));
        hold on;
    end
    title(files(ff).name(4:end-4), 'interpreter', 'none', 'fontweight', 'normal');
    drawnow;
end

% simtime=unique(traj(:,1));
% for k=1:numel(simtime)
%     figure(1); gcf; clf;
%     pos=traj(traj(:,1)==simtime(k),:);
%     plot(pos(:,3), pos(:,4), 'o');
%     axis([0 10 0 10]);
%     drawnow()
% end


function calculate_distmat()

inf1=100000;
files=dir('../PanicPackage/dump/sd_u0*100agents_*.csv');

side=10; % cm

if contains(files(1).name, '100agents')
    side=20;
end

for ff=4:size(files,1)

    fprintf('processing %s ...\n', files(ff).name);
    
    traj=csvread(['../PanicPackage/dump/', files(ff).name]);
    ids=unique(traj(:,2));
    
    k=sort(unique(traj(:,1)));
    dt=k(2)-k(1);
    
    D=zeros(numel(ids), numel(ids), numel(k));
%     tic
    for kk=1:numel(k)
        for ii=1:numel(ids)
            pos_ii=traj(traj(:,1)==k(kk) & traj(:,2)==ids(ii),3:4);     
            for jj=ii:numel(ids)
                if jj==ii
                    D(ii,jj,kk)=inf1/2;
                else
                    pos_jj=traj(traj(:,1)==k(kk) & traj(:,2)==ids(jj),3:4);
                    if ~isempty(pos_ii) && ~isempty(pos_jj)
                        pos_ii=pos_ii(1,:); % only get the first row because the same timestamp may be recorded for two positions
                        pos_jj=pos_jj(1,:); % only get the first row 
                        % if the particle is within the arena then  
                        if pos_ii(1)<side && pos_jj(1)<side
                            % calculate distance
                            try
                            D(ii,jj,kk)=norm(pos_ii-pos_jj);
                            catch
                                keyboard
                            end
                        else
                            % otherwise it's gone
                            D(ii,jj,kk)=inf1;
                        end
                    else
                        D(ii,jj,kk)=inf1;
                    end
                end
            end
        end
        D(:,:,kk)=D(:,:,kk)+D(:,:,kk)';
%         fprintf('%d, %.3f\n', kk, toc);
    end
 
    save(sprintf('../PanicPackage/dump/distmat_%s.mat', ...
        files(ff).name(1:end-4)), 'D', 'ids', 'dt', 'k');
end

function calculate_exposure_from_IA(tau)


files=dir('../PanicPackage/dump/distmat_sd_u1*100agents_*.mat');

iter=10;
side=10; % cm

if contains(files(1).name, '100agents')
    side=20;
end

for ff=1:size(files,1)

    fprintf('processing %s ...\n', files(ff).name);
    load(['../PanicPackage/dump/', files(ff).name], 'D', 'ids', 'k', 'dt');
        
    nn_dist=zeros(numel(ids),numel(k));
    for kk=1:numel(k)
        D1=D(:,:,kk);
        nn_dist(:,kk)=min(D1,[],2);
        nn_dist(nn_dist(:,kk)>side*sqrt(2),kk)=nan;
    end
    
    csvwrite(sprintf('../PanicPackage/dump/nn_dist_%s.csv', files(ff).name(9:end-4)),...
        nn_dist);
    
    max_exposure1=zeros(1,iter);
    
    for ii=1:iter
         % pick a random id as infected
        ia_id=randi(numel(ids));
        allexp=sum(exp(-tau*D),3)*dt;
        max_exposure1(ii)=max(allexp(ia_id,:));
    end
    csvwrite(sprintf('../PanicPackage/dump/mxp%.2d_%s.csv', tau, files(ff).name(9:end-4)),...
        max_exposure1);
end


function composite_plot()

% figure 
% subplot 1 and 2 eps files
% tau dependence, showing both 
% violin plots for ANND and evacuation time

lw=3; % linewidth;
fs=16; % fontsize

conditions=1:2;
tauval=[1 5 8 10 16];
% ids={'u1s1', 'u1s0', 'u0s1', 'u0s0'};
% desc={'evacuate and distance', 'evacuate only',...
%         'exit and distance', 'exit only'};
ids={'u0s1', 'u1s0'};
% ids={'u0s1_100agents', 'u1s0_100agents'};
desc={'Exit & distance', 'Evacuate'};
clr={[64,224,208]/255, [1 0 0]}; % turquoise and red


figure(1); gcf; clf;
set(gcf, 'position', [896, 1, 1019, 984]);

subplot(4,4,[1,5])
img=imread('./setup_cartoon.png');
imshow(img);
text(0, 0, '(A)', 'fontsize', fs*1.5, 'fontweight', 'bold');

subplot(4,4,[2])
% subplot(4,4,[1,2])
img=imread('../PanicPackage/dump/eps_u0s1/sd1200.jpg');
img1=imcrop(img, [650, 370, 1700, 900]);
imshow(img1);
title(desc{1}, 'fontweight', 'normal');
set(gca, 'fontsize',fs, 'units', 'normalized');
text(0, 0, '(B)', 'fontsize', fs*1.5, 'fontweight', 'bold');

subplot(4,4,[6]);
% subplot(4,4,[5,6]);
img=imread('../PanicPackage/dump/eps_u1s0/sd1083.jpg');
img1=imcrop(img, [650, 370, 1700, 900]);
imshow(img1);
title(desc{2}, 'fontweight', 'normal');
set(gca, 'fontsize',fs);
set(gca, 'fontsize',fs, 'units', 'normalized');
text(0, 0, '(C)', 'fontsize', fs*1.5, 'fontweight', 'bold');




for cc=1:numel(conditions) % for each condition e.g. u0s1
    for tt=1:numel(tauval) % for each tauvalue e.g. tau=2,3,
        mxp1=[];
        % get all the files 
        simfiles{tt}=dir(sprintf('../PanicPackage/dump/mxp%.2d*%s_*.*', tauval(tt), ids{cc}));
        for ss=1:size(simfiles{tt},1)
            mxp1=[mxp1; csvread(['../PanicPackage/dump/', simfiles{tt}(ss).name])];
        end
                
        % update dataset for a tau value for each condition
        bpdata(tt,cc,:)=mxp1(:);
    end
end

subplot(4,4,[3,4, 7,8]);
% interpolation

for cc=1:numel(conditions)
    data2=squeeze(bpdata(:,cc,:));
    
    if sum(data2<0)
        error('check!');
    end
    data20=data2; mu10=mean(data20,2);
    data2=log10(data2);
    mu1=mean(data2,2);
    std1=std(data2,[],2);
%     errorbar(tauval, mu1', std1', 'color', clr{cc});
    [hl,hp]=boundedline(tauval,mu1',std1', '-o');
    set(hl, 'color', clr{cc}, 'linewidth', lw);
    set(hp, 'facecolor', clr{cc}, 'facealpha', 0.5);

    hold on;
    
%     mu_i=interp1(tauval, mu10, tvals1, 'linear', 'extrap')
    
end
% plot(tauval, log10(900*exp(-tauval*2)), 'k--');
%xlabel('$\mathrm{m}^{-1}$', 'interpreter', 'latex');
xlabel([char(964), ' (m^{-1})']);
ylabel('Risk of exposure, E (s)');
set(gca, 'ytick', [-9:2:1]);
set(gca, 'xtick', [1:5:16]);
set(gca, 'xlim', [1 16]);
set(gca, 'yticklabel', sprintf('10^{%d}\n', get(gca, 'ytick')))
% set(gca, 'YScale', 'log');
grid on;
set(gca, 'fontsize', fs);  
set(gca, 'fontsize',fs, 'units', 'normalized');
text(0, 2, '(D)', 'fontsize', fs*1.5, 'fontweight', 'bold');





% create smaller axes in top right, and plot on it
% axes('Position',[.725 .725 .15 .15])
% box on
% dist=0:.1:3;
% lsp={'-k', ':k', '--k', '-.k'};
% for tt=[1 numel(tauval)]
%     plot(dist, exp(-tauval(tt)*dist), lsp{tt}, 'linewidth', 2);
%     hold on;
% end
% xlabel('distance (m) ');
% ylabel('probability of infecting');
% legend({'$\tau=2\ \mathrm{m}^{-1}$', '$\tau=5\ \mathrm{m}^{-1}$'}, 'interpreter', 'latex');
% grid on; 
% set(gca, 'fontsize', 12);

bpdata=[];
for cc=1:numel(conditions)
    simfiles=dir(['../PanicPackage/dump/nn_dist*', ids{cc}, '*.*']);
    nndist1=[];
    for ss=1:size(simfiles,1)
        nndist=csvread(['../PanicPackage/dump/', simfiles(ss).name]);
        nndist1(:,ss)=nanmean(nndist,2);
    end
    bpdata(:,cc)=nndist1(:);
end
subplot(4,4,[9,10,13,14]);
% bh=boxplot(bpdata);
% set(bh, 'linewidth', 2);

violins=violinplot(bpdata);
for cc=1:numel(conditions)
    violins(cc).ViolinColor=clr{cc};
    violins(cc).BoxColor='k';
    violins(cc).EdgeColor=clr{cc};
end

% distributionPlot(bpdata, 'color', clr);

xlabel('condition');
ylabel('Avg. distance to nearest neighbor (m)');
set(gca, 'xtick', 1:numel(conditions), 'xticklabel', desc)
% xtickangle(30);
grid on;
set(gca, 'fontsize', fs);  
set(gca, 'fontsize',fs, 'units', 'normalized');
text(0, 2, '(E)', 'fontsize', fs*1.5, 'fontweight', 'bold');


bpdata=[];
for cc=1:numel(conditions)
    simfiles=dir(['../PanicPackage/dump/evactime*', ids{cc}, '*.*']);
    evac_time=[];
    for ss=1:size(simfiles,1)
        evac_time=[evac_time; csvread(['../PanicPackage/dump/', simfiles(ss).name])];
    end
    bpdata(:,cc)=evac_time;
end
subplot(4,4,[11,12,15,16]);

violins=violinplot(bpdata);
for cc=1:numel(conditions)
    violins(cc).ViolinColor=clr{cc};
    violins(cc).BoxColor='k';
    violins(cc).EdgeColor=clr{cc};
end

% distributionPlot(bpdata, 'color', clr);

xlabel('condition');
ylabel('Leaving time for 23 individuals (s)');
set(gca, 'xtick', 1:numel(conditions), 'xticklabel', desc)
% xtickangle(30);
grid on;
set(gca, 'fontsize', fs);  
set(gca, 'fontsize',fs, 'units', 'normalized');
text(0, 2, '(F)', 'fontsize', fs*1.5, 'fontweight', 'bold');

