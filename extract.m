function df_mat_weighted=extract(p,T)

p_x=p(1);
p_y=p(2);
p_z=p(3);

if(p_x>64)
    p_x = 64;
end
if(p_y>64)
    p_y = 64;
end
if(p_z>50)
    p_z = 50;
end

df_vertex(1, :, :) = T(floor(p_x), floor(p_y), floor(p_z),:);
df_vertex(2, :, :) = T(floor(p_x), ceil(p_y), floor(p_z),:);
df_vertex(3, :, :) = T(ceil(p_x), floor(p_y), floor(p_z),:);
df_vertex(4, :, :) = T(ceil(p_x), ceil(p_y), floor(p_z), :);
df_vertex(5, :, :) = T(floor(p_x), floor(p_y), ceil(p_z), :);
df_vertex(6, :, :) = T(floor(p_x), ceil(p_y), ceil(p_z), :);
df_vertex(7, :, :) = T(ceil(p_x), floor(p_y), ceil(p_z),:);
df_vertex(8, :, :) = T(ceil(p_x), ceil(p_y), ceil(p_z), :);

dc = df_vertex(4, :, :) * (1 + p_x - ceil(p_x)) + df_vertex(2, :, :) * (ceil(p_x) - p_x);
dd = df_vertex(3, :, :) * (1 + p_x - ceil(p_x)) + df_vertex(1, :, :) * (ceil(p_x) - p_x);
da = dc * (1 + p_y - ceil(p_y)) + dd * (ceil(p_y) - p_y);

de = df_vertex(8, :, :) * (1 + p_x - ceil(p_x)) + df_vertex(6, :, :) * (ceil(p_x) - p_x);
df = df_vertex(7, :, :) * (1 + p_x - ceil(p_x)) + df_vertex(5, :, :) * (ceil(p_x) - p_x);
db = de * (1 + p_y - ceil(p_y)) + df * (ceil(p_y) - p_y);


df_mat_weighted(:, :) = db * (1 + p_z - ceil(p_z)) + da * (ceil(p_z) - p_z);
