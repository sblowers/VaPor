function [Tt, T_Art, T_Vein] = extractTransientResults(Timestep)

    T_TransientStore = evalin('base', 'T_TransientStore');
    DomTot = evalin('base', 'DomTot');
    Vessel1 = evalin('base', 'Vessel1');
    
    Div1 = numel(DomTot(DomTot));
    Div2 = size(Vessel1,1);
    
    Tt = zeros(size(DomTot));
    Tt(DomTot) = T_TransientStore(1:Div1,Timestep);
    Tt(~DomTot) = NaN;
    
    T_Art = T_TransientStore(Div1+1:Div2,Timestep);
    T_Vein = T_TransientStore(Div2+1:end,Timestep);

end