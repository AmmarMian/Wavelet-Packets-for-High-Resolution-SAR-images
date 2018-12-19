function Affiche_map(x,y,X, ti)

% Initialize coefficient
% Plot
% X = 20*log10(abs(peaks));
Xmax = max(X(:));
Xmin = min(X(:));

if Xmax== Xmin,
    Xmin = Xmin*(1-eps);
    Xmax = Xmax*(1+eps);
end

S.fh = figure('units','pixels',...
              'position',[644 10 794 564], 'resize', 'on');

S.axesH = axes('Units'     , 'Normalized'  , 'Position', [0.1095 0.2890 0.7921 0.6666] );
S.sliderH = uicontrol('Style'     , 'Slider'    , ...
    'Units'     , 'Normalized'    , ...
    'Position', [0.1863 0.1453 0.6372 0.0354], 'Min'       , Xmin        , ...
    'Max'       , Xmax  , 'SliderStep',[0.001 0.01], ...
    'Value'     , Xmax);
S.sliderL = uicontrol('Style'     , 'Slider'    , ...
    'Units'     , 'Normalized'    , ...
    'Position', [0.1863 0.0762 0.6372 0.0354], ...
    'Min'       ,  Xmin, ...
    'Max'       , Xmax    , ...
    'SliderStep',[0.001 0.01], ...
    'Value'     , Xmin);
S.X = X;
S.fig = imagesc(x,y,X, [Xmin Xmax]);
xlabel('y (m)')
ylabel('x (m)')
set(gca,'Ydir','normal')
colorbar
title('');
Pfa = sum(X(:) >= Xmin) / numel(X);
title([ti, sprintf(' : Dyn : [%1.1f %1.1f] dB   Pfa = %e', Xmin, Xmax, Pfa)]);
S.ti = ti;
set([S.sliderH,S.sliderL],'call',{@ed_call,S}); 

end

function [] = ed_call(varargin)
% Callback for the edit box and slider.
[h,S] = varargin{[1,3]};  % Get calling handle and structure.

switch h  % Who called?
    case S.sliderH
        Xmax = get(S.sliderH,'value');  % Get the slider's info.
        Xmin = get(S.sliderL,'value');  % Get the slider's info.
        if Xmax <= Xmin
            Xmax = Xmin + 0.00005*(get(S.sliderL, 'Max') - get(S.sliderL, 'Min'));
            set(S.sliderH,'value', Xmax);
        end
    case S.sliderL
        Xmax = get(S.sliderH,'value');  % Get the slider's info.
        Xmin = get(S.sliderL,'value');  % Get the slider's info.
        if Xmin >= Xmax
            Xmin = Xmax - 0.00005*(get(S.sliderL, 'Max') - get(S.sliderL, 'Min'));
            set(S.sliderL,'value', Xmin);
        end
    otherwise
        % Do nothing, or whatever.
end
set(S.axesH, 'Clim', [Xmin, Xmax])
Pfa = sum(S.X(:) >= Xmin) / numel(S.X);
title([S.ti, sprintf(' : Dyn : [%1.1f %1.1f] dB   Pfa = %e', Xmin, Xmax, Pfa)]);
end
