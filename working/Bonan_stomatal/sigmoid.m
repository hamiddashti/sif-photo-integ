clear
clc

% x=[0:1:10]
% y=sig(x)
% plot(x,y)
x = 1:1:10;
a_max = 10;
c=0.1:0.1:1;
c_50=0.75;

% hold on
% p=[0.5,1,2];
% cmap = colormap(autumn(length(p)))
% for i=1:length(p)
% y=sig(a_max,c,c_50,p(i));
% plot(c,y,color=cmap(i,:))
% end
% y_max=sig(a_max,c,c_50,1);
% plot(c,y_max,color="green")
% hold off

rwc = 0.1:0.01:1
A_0 = repelem(30,length(rwc));
ss = 0.1:0.5:3;

rwc_t = repelem(0.9,length(rwc));
rwc_0 = repelem(0.4,length(rwc));
cmap = jet(length(ss));

Legend=cell(length(ss),1)
figure(1)
hold on
for i = 1:length(ss)
    s= repelem(ss(i),length(rwc));
    y=arrayfun(@water_stress,rwc, A_0,s,rwc_t,rwc_0)
    
    plot(rwc,y,'*',color=cmap(i,:))
    Legend{i}=strcat('Slope:', num2str(s(i)));
    set(gca,'xDir','reverse')
    ylim([-5,35])

end
legend(Legend)
xlabel("RWC",FontSize=12)
ylabel("Photosynthesis",FontSize=12)
saveas(gcf,"./Figures/piecewise.png")
% y=piece(A_0,s,rwc,rwc_t,rwc_0)
function [y]=sig(a_max,c,c_50,p)
    %Based on "Analyses of crp salt tolerance data"
    y = a_max./(1+((1-c)./c_50).^p);
end

function [y]= sig2(x)
    y = 1./(1+1./exp(-x))
end

function y=water_stress(rwc,A_0,s,rwc_t,rwc_0)
    if rwc >rwc_t
        y=A_0;
    elseif (rwc<=rwc_t)&&(rwc>rwc_0)
        y = A_0-A_0.*s.*(rwc_t-rwc)
    elseif (rwc <= rwc_0)
        y=0;
    end

end
