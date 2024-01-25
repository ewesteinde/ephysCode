            function subplotsqueezeV(fig, nF)
                hAx = findobj(fig, 'type', 'axes');
                for h = 1:length(hAx)
                vCurrPos = get(hAx(h), 'position'); % current position
                set(hAx(h), 'position', (vCurrPos.*[1 1 1 nF])-[0 vCurrPos(4)*(nF-1)/2 0 0]);
                end
            end