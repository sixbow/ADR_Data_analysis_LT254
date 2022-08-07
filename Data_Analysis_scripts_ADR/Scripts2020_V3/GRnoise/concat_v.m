function vec_out = concat_v(iter,object,dim)
vec_out = object{1,1};
if dim==1
j = 1;
for i=iter
   vec_out = [vec_out object{i,j}];
end  
elseif dim==2
i = 1;
for j=iter
   vec_out = [vec_out object{i,j}];
end  
else
    error('Paniek Paniek!! Hier kan ik niets mee!')
end
end

