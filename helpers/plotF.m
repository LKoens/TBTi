function plotF(f,S,Phi)

    numS = size(S,1);
    numPhi = size(S,2);

    figure
    tiledlayout(6,3)
    for i = 1 : 6
        F = reshape(f(:,i), 3, []);
        Fx = reshape(F(1,:), numS, numPhi);
        Fy = reshape(F(2,:), numS, numPhi);
        Fz = reshape(F(3,:), numS, numPhi);

        nexttile()
        surf(S,Phi,Fx)
        view(0,90)
        shading flat
        xlabel('$s$')
        ylabel('$\phi$')
        c = colorbar;
        c.TickLabelInterpreter = 'latex';
        colormap(viridis)

        nexttile()
        surf(S,Phi,Fy)
        view(0,90)
        shading flat
        xlabel('$s$')
        ylabel('$\phi$')
        c = colorbar;
        c.TickLabelInterpreter = 'latex';
        colormap(viridis)

        nexttile()
        surf(S,Phi,Fz)
        view(0,90)
        shading flat
        xlabel('$s$')
        ylabel('$\phi$')
        c = colorbar;
        c.TickLabelInterpreter = 'latex';
        colormap(viridis)
    end
end