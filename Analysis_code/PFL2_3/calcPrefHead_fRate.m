function prefHead = calcPrefHead_fRate(bData, tData)

    total_mov_mm = abs(bData.vel_for + abs(bData.vel_side) + abs(deg2rad(bData.vel_yaw)*4.5));
    no0vel_idx = find(total_mov_mm > 0);
    angle = bData.angle(no0vel_idx); 

    activity = tData.fRate_sec;
    activity = activity(no0vel_idx);

    edges_angle = [-180:60:180];
    [angleBins, centers_angle, ~] = binData(activity(~isnan(activity)), angle(~isnan(angle)), edges_angle);
    prefHead_smooth = smoothdata(angleBins,'gaussian',1);
    figure(11);clf;
    plot(centers_angle,angleBins)
    hold on
    plot(centers_angle,prefHead_smooth)
    %if sum(isnan(prefHead_smooth))
       prefHead = centers_angle(prefHead_smooth == max(prefHead_smooth));
%     else
%         bump_params = fit_sinusoid(prefHead_smooth, [0,2*pi], 1);
%         prefHead = wrapTo180((rad2deg(wrapToPi(bump_params.pos_rad + deg2rad(15)))) - 180);
%     end
end