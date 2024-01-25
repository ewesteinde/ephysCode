activity_values = [];
x_values = [];
y_values = [];
angle = [];
for t = 2:length(processed_trialData)
        
    activity_values = cat(1,activity_values, processed_trialData{t}.fRate_sec(0.5*1000:end-(0.5*1000)));
    x_values = cat(1, x_values, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
    y_values = cat(1, y_values, processed_behaviourData{t}.vel_for(0.5*1000:end-(0.5*1000)));
    angle = cat(1, angle, processed_behaviourData{t}.angle(0.5*1000:end-(0.5*1000)));
end


    %activity_values = normalize(activity_values);

    x_edges = [-10:1:10];
    y_edges = [-1:1:8];

%activity_values = normalize(activity_values);

[N, heatmap_array, x_centers, y_centers] = create_activity_heatmap(x_values, y_values, activity_values, x_edges, y_edges);

figure(1);clf;

    g(1) = subplot(4,1,1);
    plot(angle)
    ylabel('angle')
    title('111120 Cell 2')

    g(2) = subplot(4,1,2);
    plot(x_values)
    ylabel('side mm/s')

    
    g(3) = subplot(4,1,3);
    plot(y_values)
    ylabel('for mm/s')
    
    g(4) = subplot(4,1,4);
    plot(activity_values)
    ylabel('spikes/sec')
    
    linkaxes(g,'x')
    
    

figure(2); clf;
            imagesc([-10 10], [],flip(heatmap_array))
            colorbar
            title('111120 Cell 2')



