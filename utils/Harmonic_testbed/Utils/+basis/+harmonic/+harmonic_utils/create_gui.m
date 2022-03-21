function create_gui(S1,S2,X1,X2,B1,B2)
    p = uipanel('Title','Controls',...
             'Position',[0 0 1 1]);
         
%     uicontrol(p,'style','text','string','IniSize:','Units','normalized',...
%         'Position',[0 0.9 0.5 0.1]);
%     uicontrol(p,'style','Edit','string','4', 'Units','normalized',...
%          'Tag','inisize','Position',[0.6 0.9 0.25 0.1]);
%     
%     uicontrol(p,'style','text','string','FinalSize:','Units','normalized',...
%         'Position',[0 0.75 0.5 0.1]);
%     uicontrol(p,'style','Edit','string','40', 'Units','normalized',...
%          'Tag','largesize','Position',[0.6 0.75 0.25 0.1]);
% 
%     uicontrol(p,'style','text','string','Step:','Units','normalized',...
%         'Position',[0 0.6 0.5 0.1]);
%     uicontrol(p,'style','Edit','string','5', 'Units','normalized',...
%          'Tag','step','Position',[0.6 0.6 0.25 0.1]);
%      
%     uicontrol(p,'style','text','string','Reg weight:','Units','normalized',...
%         'Position',[0 0.4 0.5 0.15]);
%     uicontrol(p,'style','Edit','string','0.1', 'Units','normalized',...
%          'Tag','reg_weight','Position',[0.6 0.4 0.3 0.15]);
%      
%     uicontrol(p,'style','text','string','Ini type:','Units','normalized',...
%         'Position',[0 0.2 0.5 0.15]);
%     uicontrol(p,'style','popupmenu','string',{'nearest neigh','landmarks'},...
%         'Units','normalized', 'Value', 1, 'Tag','ini_type','Position',[0.5 0.2 0.5 0.15]);
     
    c = uicontrol(p,'Units','normalized','Position',[0 0.9 0.9 0.1]);
    c.String = 'Clear Landmarks';
    c.Callback = @clear_landmarks;
        
    c = uicontrol(p,'Units','normalized','Position',[0 0.7 0.9 0.1]);
    c.String = 'Harmonic Matching';
    c.Callback = {@harmonic_interact,S1,S2,X1,X2,B1,B2};
    
%     c = uicontrol(p,'Units','normalized','Position',[0 0.6 0.9 0.1]);
%     c.String = 'Power Harmonic Matching';
%     c.Callback = {@opt_harmonic_interact,S1,S2,X1,X2,B1,B2};
    
    c = uicontrol(p,'Units','normalized','Position',[0 0.6 0.9 0.1]);
    c.String = 'Central FE Optimization';
    c.Callback = {@central_FE_interact,S1,S2,X1,X2,B1,B2};

    c = uicontrol(p,'Units','normalized','Position',[0 0.5 0.9 0.1]);
    c.String = 'Sqrt3 Intrinsic Refinement';
    c.Callback = {@double_intrinsic_sqrt3,S1,S2,X1,X2,B1,B2};
    
    c = uicontrol(p,'Units','normalized','Position',[0 0.4 0.9 0.1]);
    c.String = 'Sqrt3 Power Refinement';
    c.Callback = {@double_intrinsic_sqrt3_power,S1,S2,X1,X2,B1,B2};
    
    
%     c = uicontrol(p,'Units','normalized','Position',[0 0.5 0.9 0.1]);
%     c.String = 'Double Power Harmonic Matching';
%     c.Callback = {@double_opt_harmonic_interact,S1,S2,X1,X2,B1,B2};
    
%     c = uicontrol(p,'Units','normalized','Position',[0 0.4 0.9 0.1]);
%     c.String = 'Iterative Power Harmonic Matching';
%     c.Callback = {@iterative_opt_harmonic_interact,S1,S2,X1,X2,B1,B2};
    
%     c = uicontrol(p,'Units','normalized','Position',[0 0.3 0.9 0.1]);
%     c.String = 'Double Iterative Power Harmonic Matching';
%     c.Callback = {@double_iterative_opt_harmonic_interact,S1,S2,X1,X2,B1,B2};
    
    c = uicontrol(p,'Units','normalized','Position',[0 0.15 0.9 0.1]);
    c.String = 'Control Landmark Double';
    c.Callback = {@double_control_landmarks,S1,S2,X1,X2,B1,B2};
    
    c = uicontrol(p,'Units','normalized','Position',[0 0 0.9 0.1]);
    c.String = 'Refinement Control Landmark';
    c.Callback = {@double_control_refinement,S1,S2,X1,X2,B1,B2};
    
%     c = uicontrol(p,'Units','normalized','Position',[0 0 0.9 0.1]);
%     c.String = 'Fourier Harmonic Matching';
%     c.Callback = {@fourier_harmonic_interact,S1,S2,X1,X2,B1,B2};
   
    
    
    
end