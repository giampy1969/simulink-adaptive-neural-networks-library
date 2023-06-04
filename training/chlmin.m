
tic

% open scheme
load_system('chlnet');

load choles_all
t=t(:,1:130);p=p(:,1:130);
mindata=min([t' p']);maxdata=max([t' p']);

% nonlinear network gain and switch
set_param('chlnet/Switch','sw','1');

% nonlinear network initial conditions
nwrap=str2num(get_param('chlnet/GRBF','nwrap'));
set_param('chlnet/GRBF','x0','W(end,:)');
W=zeros(((21+2)*100+2)*3,1)';

for dgen=[0.6 0.5 0.4 0.3 0.2],
    
    str=[ '[' num2str(dgen) '  ' num2str(dgen) '  1]' ]
    set_param('chlnet/GRBF','e2',str);
    
    for i=1:50000,sim('chlnet');end
    
    pause(0.5);drawnow;

    W=W(end,:);
    str=['save chl_W' num2str(fix(dgen*10)) ' W mindata maxdata -v4' ];
    eval(str)
    
end

pause(0.5);
load choles_all
t=t(:,131:end);p=p(:,131:end);

for d=6:-1:2,
    eval(['load chl_W' num2str(d) ]);

    pause(0.5);drawnow;

    sim('chlnet');

    for i=1:3,
        per_rms=100*sqrt((est(:,i)-nom(:,i)).^2)./(ones(length(est),1)*(max(nom(:,i))-min(nom(:,i))));
        mean_per_rms(i)=mean(per_rms);
    end;

    sum(mean_per_rms)

end

toc/3600
    
