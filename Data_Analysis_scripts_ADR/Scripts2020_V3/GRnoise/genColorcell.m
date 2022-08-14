function Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors)
%Colorcell = genColorcell(T_iter,Tindex,Tcolors,P_iter,Pindex,Pcolors)
if length(T_iter)>1
Colorcell = [Tcolors(Tindex,:),{'black'},{'black'},Tcolors(Tindex,:),Tcolors(Tindex,:)];
elseif length(P_iter)>1
Colorcell = [Pcolors(Pindex,:),{'black'},{'black'},Pcolors(Pindex,:),Pcolors(Pindex,:)];
else
Colorcell = [Tcolors(Tindex,:),{'black'},{'black'},Tcolors(Tindex,:),Tcolors(Tindex,:)];
end
end

