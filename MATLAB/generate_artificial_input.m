% Compute artificial alpha input (10 Hz sine wave) 
% First part:  ~10 minutes of activity with short alpha burst
% Second part: ~10 minutes of activity with slow power
% modulation (0.01 - 0.03 Hz)


function generate_artificial_input()

    %{
    Generate first part: 10 minutes with short alpha burst
    %}

    % Generate one alpha cycle @10000 Hz
    sp          =   linspace(0,2*pi,1001);
    onecycle    =   sin(sp);
    onecycle    =   zscore(onecycle(1:end-1));

    % Generate four seconds of alpha activity
    tinsec      =   4;
    inp         =   [];
    for i = 1:10*tinsec
        inp     =   [inp onecycle];
    end
    
    % Generate a slow modulator of the power of the 10 Hz sine wave
    tslen       =   5000;
    mod1        =   linspace(0.005,0.4,tslen);
    mod2        =   linspace(0.4,0.001,tslen);
    hcyc0       =   linspace(pi,2*pi,tslen+1);
    hcyc1       =   linspace(pi,2*pi,tslen+1);
    hcyc2       =   linspace(2*pi,3*pi,tslen+1);
    hcyc0       =   ones(1,length(hcyc0)-1)*0.02;
    hcyc1       =   zscore(cos(hcyc1(1:end-1))) * std(mod1) + mean(mod1);
    hcyc2       =   zscore(cos(hcyc2(1:end-1))) * std(mod2) + mean(mod2);
    osc_modulator = [hcyc0 hcyc0 hcyc0 hcyc1 hcyc2 hcyc0 hcyc0 hcyc0];

    % Compute alpha burst
    ext_input   =   zscore(inp(10000:30000).*osc_modulator(10000:30000));

    % Generate 10 minutes of 10 Hz sine wave with the power burst in the middle
    ext_input   =   [repmat(ext_input(1:1000),1,5*60*10) ext_input repmat(ext_input(1:1000),1,5*60*10)];
    cut_sig     =   downsample(repmat(ext_input(1:1000),1,10*60*10),10);
    
    % downsample to 1000 Hz
    ext_input   =   downsample(ext_input,10);
    
    
    
    %{
    Generate second part: 10 minutes with slow power modulation
    %}

    % Generate 310 TRs of a 10 Hz sine wave (1 TR = 1.94 sec)
    tinsec      =   1.94*310;
    sp          =   linspace(0,2*pi*tinsec*10,tinsec*1000);
    alphawave   =   sin(sp);

    % Generate power modulators at different frequencies
    sp2_p01     =   linspace(0,2*pi,1000/0.01);
    sp2_p02     =   linspace(0,2*pi,1000/0.02);
    sp2_p03     =   linspace(0,2*pi,1000/0.03);
    sp2_p015    =   linspace(0,2*pi,1000/0.015);
    sp2_p025    =   linspace(0,2*pi,1000/0.025);
    modulator   =   cos([sp2_p01 sp2_p02 sp2_p03 sp2_p03 sp2_p01 sp2_p025 sp2_p025 sp2_p015 sp2_p01 sp2_p03 sp2_p03]);
    modulator   =   modulator + -1*min(modulator);
    inp_seg2    =   zscore(alphawave.*modulator(1:length(alphawave)));

    % Concat first and second part
    ext_input   =   [ext_input inp_seg2];
    
    % Fill up edges to get 668 TRs of activity (to match length of
    % empirical data)
    missing     =   668*1.94*1000 - length(ext_input);
    ext_input   =   [cut_sig(1:round(missing/2)) ext_input cut_sig(1:round(missing/2)+1000)];

    % Copy for all 68 nodes of model
    output      =   repmat(ext_input',1,68);

    dlmwrite('Arificial_input.txt',output,'delimiter',' ','precision','%.4f');
end

