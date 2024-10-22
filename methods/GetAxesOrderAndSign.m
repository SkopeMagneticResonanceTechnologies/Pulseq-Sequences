function [axesOrder, axesSign, readDir_SCT, phaseDir_SCT, sliceDir_SCT] = GetAxesOrderAndSign(sliceOrientation,phaseEncDir, doFlipXAxis)
% SIEMENS physical coordinate system: 
%          X: (Your) Left -> (Your) Right   |   Y:        Down -> Up     | Z:  Rear -> Front (Towards you)  
%
% Vendor independent patient coordinte system: SAG/COR/TRA for HFS 
%          S: (Your) Left -> (Your) Right   |   C:        Up -> Down     | T:  Front -> Rear (Away from you)      

% (c) 2024 Skope Magnetic Resonance Technologies AG

        switch sliceOrientation
                case SliceOrientation.TRA
                    switch phaseEncDir
                        case PhaseEncodingDirection.AP
                            axesOrder = {'x','y','z'};
                            axesSign = [1,1,1];
                            readDir_SCT = [1,0,0];
                            phaseDir_SCT = [0,1,0];
                            sliceDir_SCT = [0,0,1];
                        case PhaseEncodingDirection.PA
                            axesOrder = {'x','y','z'};
                            axesSign = [1,-1,1];
                            readDir_SCT = [1,0,0];
                            phaseDir_SCT = [0,-1,0];
                            sliceDir_SCT = [0,0,1];
                        case PhaseEncodingDirection.RL
                            axesOrder = {'y','x','z'};
                            axesSign = [1,-1,1];
                            readDir_SCT = [0,-1,0];
                            phaseDir_SCT = [1,0,0];
                            sliceDir_SCT = [0,0,1];
                        case PhaseEncodingDirection.LR
                            axesOrder = {'y','x','z'};
                            axesSign = [1,1,1];
                            readDir_SCT = [0,-1,0];
                            phaseDir_SCT = [-1,0,0];
                            sliceDir_SCT = [0,0,1];
                        otherwise
                            error('Wrong phase encoding direction.');
                    end
                case SliceOrientation.SAG   
                    switch phaseEncDir
                        case PhaseEncodingDirection.AP
                            axesOrder = {'z','y','x'};
                            axesSign = [1,1,1];
                            readDir_SCT = [0,0,-1];
                            phaseDir_SCT = [0,1,0];
                            sliceDir_SCT = [-1,0,0];
                        case PhaseEncodingDirection.PA
                            axesOrder = {'z','y','x'};
                            axesSign = [1,-1,1];
                            readDir_SCT = [0,0,-1];
                            phaseDir_SCT = [0,-1,0];
                            sliceDir_SCT = [-1,0,0];
                        case PhaseEncodingDirection.HF
                            axesOrder = {'y','z','x'};
                            axesSign = [1,-1,1];
                            readDir_SCT = [0,-1,0];
                            phaseDir_SCT = [0,0,-1];
                            sliceDir_SCT = [-1,0,0];
                        case PhaseEncodingDirection.FH
                            axesOrder = {'y','z','x'};
                            axesSign = [1,1,1];
                            readDir_SCT = [0,-1,0];
                            phaseDir_SCT = [0,0,1];
                            sliceDir_SCT = [-1,0,0];
                        otherwise
                            error('Wrong phase encoding direction.');
                    end
                case SliceOrientation.COR
                    switch phaseEncDir
                        case PhaseEncodingDirection.HF
                            axesOrder = {'x','z','y'};
                            axesSign = [1,-1,1];
                            readDir_SCT = [1,0,0];
                            phaseDir_SCT = [0,0,-1];
                            sliceDir_SCT = [0,1,0];
                        case PhaseEncodingDirection.FH
                            axesOrder = {'x','z','y'};
                            axesSign = [1,1,1];
                            readDir_SCT = [1,0,0];
                            phaseDir_SCT = [0,0,1];
                            sliceDir_SCT = [0,1,0];
                        case PhaseEncodingDirection.RL
                            axesOrder = {'z','x','y'};
                            axesSign = [1,-1,1];
                            readDir_SCT = [0,0,-1];
                            phaseDir_SCT = [1,0,0];
                            sliceDir_SCT = [0,1,0];
                        case PhaseEncodingDirection.LR
                            axesOrder = {'z','x','y'};
                            axesSign = [1,1,1];
                            readDir_SCT = [0,0,-1];
                            phaseDir_SCT = [-1,0,0];
                            sliceDir_SCT = [0,1,0];
                        otherwise
                            error('Wrong phase encoding direction.');
                    end                    
                otherwise
                    error('Wrong slice orientation.');
    end

    if doFlipXAxis % To correct for flipped X-axis in Pulseq
        ind = find(cellfun(@(x)strcmp(x,'x'),axesOrder));
        axesSign(ind) = -axesSign(ind);
    end

end