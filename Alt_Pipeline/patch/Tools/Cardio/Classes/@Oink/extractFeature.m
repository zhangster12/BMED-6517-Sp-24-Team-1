function [ftVector, obj] = extractFeature(obj, feature, level, varargin)

% -------------------------------------------------------------------------
% This function extracts the specified feature from the specified level. If
% the data for the feature is not available, it is first computed.
%
% Arguments (required)
% - feature     Feature     Feature for extraction
% - level       Level       Level for which to extract feature
%
% Arguments (optional)
% - subject                 Subject for which to extract feature (default 1)
% -------------------------------------------------------------------------

% Parse optional input arguments
if ~isempty(varargin)
    for arg = 1:length(varargin)
        switch varargin{arg}
            case 'subject'; subject = varargin{arg + 1};
        end
    end
end

% Set defaults for optional arguments
if ~exist('subject', 'var'); subject = 1; end

% Determine the parent for this level
parent = obj.getParent(level);

% If the level and its parent do not exist, return
if ~isfield(obj.dataset{subject}, string(parent)) && ~isfield(obj.dataset{subject}, string(level))
    disp("-> Error in Oink.extractFeature(): Level." + string(level) + "does not exist");
    ftVector = nan; return
end

% Specify whether the parent for the current level exists
parentFLAG = isfield(obj.dataset{subject}, string(parent));

% Extract the feature at the specified level
switch feature
    
    case Feature.scgPEP
                    
        % SCG-derived pre-ejection period
        [obj, ftVector] = getFeature(obj, 'scgRAO', 'scgAortic', 'AO', []);

    case Feature.scgLVET

        % SCG-derived LVET
        [obj, rao] = getFeature(obj, 'scgRAO', 'scgAortic', 'AO', []);
        [obj, rac] = getFeature(obj, 'scgRAC', 'scgAortic', 'AC', []);
        ftVector = rac - rao;

    case Feature.truePEP

        % Aortic pressure-derived PEP
        [obj, ftVector] = getFeature(obj, 'trueRAO', 'trueAortic', 'AO', []);

    case Feature.trueLVET

        % Aortic pressure-derived LVET
        [obj, rao] = getFeature(obj, 'trueRAO', 'trueAortic', 'AO', []);
        [obj, rac] = getFeature(obj, 'trueRAC', 'trueAortic', 'AC', []);
        ftVector = rac - rao;

    case Feature.scgPAT

        % SCG- and PPG-derived pulse arrival time
        [obj, ftVector] = getFeature(obj, 'scgPAT', 'scgPTT', [], []);
        
    case Feature.scgPATCalibrated
        
        % SCG- and PPG-derived pulse arrival time (calibrated)
        [obj, ftVector] = getFeature(obj, 'scgPATCalibrated', 'calibratePTT', [], []);    

    case Feature.scgPTT

        % SCG- and PPG-derived pulse transit time
        [obj, ftVector] = getFeature(obj, 'scgPTT', 'scgPTT', [], []);
        
    case Feature.scgPTTCalibrated
        
        % SCG- and PPG-derived pulse transit time (calibrated)
        [obj, ftVector] = getFeature(obj, 'scgPTTCalibrated', 'calibratePTT', [], []);

    case Feature.truePAT

        % Aortic and femoral pressure-derived pulse arrival time
        [obj, ftVector] = getFeature(obj, 'truePAT', 'truePTT', [], []);
        
    case Feature.truePATCalibrated
        
        % Aortic and femoral pressure-derived pulse arrival time (calibrated)
        [obj, ftVector] = getFeature(obj, 'truePATCalibrated', 'calibratePTT', [], []);

    case Feature.truePTT

        % Aortic and femoral pressure-derived pulse transit time
        [obj, ftVector] = getFeature(obj, 'truePTT', 'truePTT', [], []);
        
    case Feature.truePTTCalibrated
        
        % Aortic and femoral pressure-derived pulse transit time (calibrated)
        [obj, ftVector] = getFeature(obj, 'truePTTCalibrated', 'calibratePTT', [], []);

    case Feature.hr

        % Heart rate (Biopac ECG)
        if parentFLAG && isfield(obj.dataset{subject}.(string(parent)), 'HR_biopac')
            ftVector = obj.getSubinterval('HR_biopac', level, 'subjects', subject);
        else; ftVector = obj.dataset{subject}.(string(level)).HR_biopac;
        end

    case Feature.hrvDifference

        % Heart rate variability (difference method)
        [obj, ftVector] = getFeature(obj, 'difference', 'ecgHRV', 'difference', 'hrv');
        
    case Feature.hrvSpectral

        % Heart rate variability (spectral method)
        [obj, lf] = getFeature(obj, 'spectralLF', 'ecgHRV', 'spectral', 'hrv');
        [obj, hf] = getFeature(obj, 'spectralHF', 'ecgHRV', 'spectral', 'hrv');
        ftVector = lf(:)./hf(:);
        
    case Feature.hrvSpectralLF
        
        % Heart rate variability (spectral, low-frequency)
        [obj, ftVector] = getFeature(obj, 'spectralLF', 'ecgHRV', 'spectral', 'hrv');
        
    case Feature.hrvSpectralHF
        
        % Heart rate variability (spectral, high-frequency)
        [obj, ftVector] = getFeature(obj, 'spectralHF', 'ecgHRV', 'spectral', 'hrv');
        
    case Feature.hrvPoincare

        % Heart rate variability (Poincare method)
        [obj, e1] = getFeature(obj, 'poincareEig1', 'ecgHRV', 'poincare', 'hrv');
        [obj, e2] = getFeature(obj, 'poincareEig2', 'ecgHRV', 'poincare', 'hrv');
        ftVector = e1(:)./e2(:);
        
    case Feature.hrvPoincareEig1
        
        % Heart rate variability (Poincare, first eigenvalue)
        [obj, ftVector] = getFeature(obj, 'poincareEig1', 'ecgHRV', 'poincare', 'hrv');
        
    case Feature.hrvPoincareEig2
        
        % Heart rate variability (Poincare, second eigenvalue)
        [obj, ftVector] = getFeature(obj, 'poincareEig2', 'ecgHRV', 'poincare', 'hrv');
        
    case Feature.scgAmplitudeApex

        % Amplitude of SCG at apex
        [obj, ftVector] = getFeature(obj, 'apexSCG_z', 'amplitudes', Mode.apexSCG_z, 'amplitude');
        
    case Feature.scgAmplitudeSternum

        % Amplitude of SCG at sternum
        [obj, ftVector] = getFeature(obj, 'sternumSCG_z', 'amplitudes', Mode.sternumSCG_z, 'amplitude');

    case Feature.ppgAmplitudeProximal

        % Amplitude of proximal PPG
        [obj, ftVector] = getFeature(obj, 'apexPPG', 'amplitudes', Mode.apexPPG, 'amplitude');

    case Feature.ppgAmplitudeDistal

        % Amplitude of distal PPG
        [obj, ftVector] = getFeature(obj, 'femoralPPG', 'amplitudes', Mode.femoralPPG, 'amplitude');

    case Feature.scgPEPOverLVET

        % SCG-derived PEP/LVET ratio
        [obj, pep] = getFeature(obj, 'scgRAO', 'scgAortic', 'AO', []);
        [obj, rac] = getFeature(obj, 'scgRAC', 'scgAortic', 'AC', []);
        ftVector = pep./(rac - pep);

    case Feature.truePEPOverLVET

        % Aortic pressure-derived PEP/LVET ratio
        [obj, pep] = getFeature(obj, 'trueRAO', 'trueAortic', 'AO', []);
        [obj, rac] = getFeature(obj, 'trueRAC', 'trueAortic', 'AC', []);
        ftVector = pep./(rac - pep);

    case Feature.proximalSPO2

        % SPO2 from proximal sensor
        if parentFLAG && isfield(obj.dataset{subject}.(string(parent)), 'apexSPO2')
            ftVector = mean(obj.getSubinterval('apexSPO2', level, 'subjects', subject));
        else; ftVector = mean(obj.dataset{subject}.(string(level)).apexSPO2);
        end

    case Feature.distalSPO2

        % SPO2 from distal sensor
        if parentFLAG && isfield(obj.dataset{subject}.(string(parent)), 'femoralSPO2')
            ftVector = mean(obj.getSubinterval('femoralSPO2', level, 'subjects', subject));
        else; ftVector = mean(obj.dataset{subject}.(string(level)).femoralSPO2);
        end

    case Feature.scgSQI

        % If the SQI is requested but not available, throw an
        % error since the SQI takes a while to generate
        if ~isfield(obj.dataset{subject}.(string(parent)), 'sqi')
            disp("-> Error in Oink.trainingData(): No SCG SQI data available"); return
        elseif ~isfield(obj.dataset{subject}.(string(parent)), string(Mode.sternumSCG_z))
            disp("-> Error in Oink.trainingData(): No SCG SQI data available"); return
        else
            ftVector = obj.getSubinterval('sternumSCG_z', level, 'subjects', subject, 'under', "sqi");
        end

    case Feature.ppgSQI

        % If the SQI is requested but not available, throw an
        % error since the SQI takes a while to generate
        if ~isfield(obj.dataset{subject}.(string(parent)), 'sqi')
            disp("-> Error in Oink.trainingData(): No PPG SQI data available"); return
        elseif ~isfield(obj.dataset{subject}.(string(parent)), string(Mode.femoralPPG))
            disp("-> Error in Oink.trainingData(): No PPG SQI data available"); return
        else
            ftVector = obj.getSubinterval('femoralPPG', level, 'subjects', subject, 'under', "sqi");
        end

    case Feature.meanAorticPressure

        % Mean aortic pressure from catheter
        [obj, ftVector] = getFeature(obj, 'aorticMAP', 'trueMAP', 'aortic', []);

    case Feature.meanFemoralPressure

        % Mean femoral pressure from catheter
        [obj, ftVector] = getFeature(obj, 'femoralMAP', 'trueMAP', 'femoral', []);

    case Feature.meanWedgePressure

        % Mean wedge pressure from catheter
        [obj, ftVector] = getFeature(obj, 'wedgeMAP', 'trueMAP', 'wedge', []);

    case Feature.meanRAPressure

        % Mean right atrial pressure from catheter
        [obj, ftVector] = getFeature(obj, 'rightAtrialMAP', 'trueMAP', 'ra', []);

    case Feature.aorticPressureSQI

        % If the SQI is requested but not available, throw an
        % error since the SQI takes a while to generate
        if ~isfield(obj.dataset{subject}.(string(parent)), 'sqi')
            disp("-> Error in Oink.trainingData(): No Aortic Pressure SQI data available"); return
        elseif ~isfield(obj.dataset{subject}.(string(parent)), string(Mode.aorticPressure))
            disp("-> Error in Oink.trainingData(): No Aortic Pressure SQI data available"); return
        else
            ftVector = obj.getSubinterval('aorticPressure', level, 'subjects', subject, 'under', "sqi");
        end

    case Feature.femoralPressureSQI

        % If the SQI is requested but not available, throw an
        % error since the SQI takes a while to generate
        if ~isfield(obj.dataset{subject}.(string(parent)), 'sqi')
            disp("-> Error in Oink.trainingData(): No Femoral Pressure SQI data available"); return
        elseif ~isfield(obj.dataset{subject}.(string(parent)), string(Mode.femoralPPG))
            disp("-> Error in Oink.trainingData(): No Femoral Pressure SQI data available"); return
        else
            ftVector = obj.getSubinterval('femoralPressure', level, 'subjects', subject, 'under', "sqi");
        end

    case Feature.systolicPressure

        % Systolic pressure from catheter
        [obj, ftVector] = getFeature(obj, 'systolic', 'trueMAP', 'aortic', []);

    case Feature.diastolicPressure

        % Diastolic pressure from catheter
        [obj, ftVector] = getFeature(obj, 'diastolic', 'trueMAP', 'aortic', []);

    case Feature.pulsePressure

        % Pulse pressure from catheter
        [obj, ftVector] = getFeature(obj, 'aorticPressure', 'amplitudes', Mode.aorticPressure, 'amplitude');
        
    case Feature.trueIHAT
        
        % iHAT: PAT/RR interval from catheter
        [obj, pat] = getFeature(obj, 'truePAT', 'truePTT', [], []);
        if parentFLAG && isfield(obj.dataset{subject}.(string(parent)), 'HR_biopac')
            hr = obj.getSubinterval('HR_biopac', level, 'subjects', subject);
        else; hr = obj.dataset{subject}.(string(level)).HR_biopac;
        end; ftVector = pat(:)./hr(:);
        
    case Feature.ppgIHAT
        
        % iHAT: PAT/RR interval from PPG
        [obj, pat] = getFeature(obj, 'scgPAT', 'scgPTT', [], []);
        if parentFLAG && isfield(obj.dataset{subject}.(string(parent)), 'HR_biopac')
            hr = obj.getSubinterval('HR_biopac', level, 'subjects', subject);
        else; hr = obj.dataset{subject}.(string(level)).HR_biopac;
        end; ftVector = pat(:)./hr(:);
        
    case Feature.femoralPressurePPV
        
        % Pulse pressure variability (femoral artery catheter)
        [obj, ftVector] = getFeature(obj, 'femoralPressure', 'respvar', Mode.femoralPPG, 'respvar');
        
    case Feature.femoralPPGPPV
        
        % Pulse pressure variability (femoral artery PPG)
        [obj, ftVector] = getFeature(obj, 'femoralPPG', 'respvar', Mode.femoralPPG, 'respvar');
    
end

% Sub-function for extracting feature and performing error checking
    function [obj, feature] = getFeature(obj, fieldName, functionName, functionFlag, under)
        
        % This function extracts a feature under the following scenarios:
        % (1) If the parent exists, but the necessary data has not been
        % extracted in the parent or current level, extract it in the
        % parent and get the sub-interval for the current level data.
        % (2) If the parent exists and contains the extracted data, obtain
        % the current level as a sub-interval of the parent
        % (3) If the data is extracted in the current level, return the
        % data from the current level
        % (4) If the data is nowhere and there is no parent, extract the
        % data into the current level and return the current level data
        
        % If the data is a direct sub-field of "obj.dataset" (i.e. the
        % field is located in obj.dataset.fieldName):
        
        if isempty(under)
        
            % Scenario (1)
            if parentFLAG && ~isfield(obj.dataset{subject}.(string(parent)), fieldName) && ...
                    ~isfield(obj.dataset{subject}.(string(level)), fieldName)
                if ~isempty(functionFlag)
                    obj = obj.(functionName)('subjects', subject, 'levels', parent, functionFlag);
                else
                    obj = obj.(functionName)('subjects', subject, 'levels', parent);
                end; feature = obj.getSubinterval(fieldName, level, 'subjects', subject);

                % Scenario (2)
            elseif parentFLAG && isfield(obj.dataset{subject}.(string(parent)), fieldName)
                feature = obj.getSubinterval(fieldName, level, 'subjects', subject);

                % Scenario (3)
            elseif isfield(obj.dataset{subject}.(string(level)), fieldName)
                feature = obj.dataset{subject}.(string(level)).(fieldName);

                % Scenario (4)
            else
                if ~isempty(functionFlag)
                    obj = obj.(functionName)('subjects', subject, 'levels', level, functionFlag);
                else
                    obj = obj.(functionName)('subjects', subject, 'levels', level);
                end; feature = obj.dataset{subject}.(string(level)).(fieldName);
            end
        
        end
        
        % If the data is not a direct sub-field of "obj.dataset" (i.e. the
        % field is located in obj.dataset.under.fieldName):
        
        % The four scenarios are updated as follows:
        % (1) If the parent exists, but the "under" field does not exist
        % for the parent or current level, extract the data into the
        % parent and return the data for the current level as a
        % sub-interval of the parent.
        % (2) If the parent exists and contains the data, return the data
        % for the current level as a sub-interval of the parent.
        % (3) If the parent doesn't have the data but the current level
        % does, return the data for the current level.
        % (4) If the parent doesn't exist and the data has not been
        % extracted in the current level, extract the data for the current
        % level and return
        
        if ~isempty(under)
        
            % Scenario (1)
            if parentFLAG && ~isfield(obj.dataset{subject}.(string(parent)), under) && ...
                    ~isfield(obj.dataset{subject}.(string(level)), under)
                obj = obj.(functionName)('subjects', subject, 'levels', parent, 'modalities', functionFlag);
                feature = obj.getSubinterval(fieldName, level, 'subjects', subject, 'under', under);

                % Scenario (2)
            elseif parentFLAG && isfield(obj.dataset{subject}.(string(parent)), under) && ...
                    isfield(obj.dataset{subject}.(string(parent)).(under), fieldName)
                feature = obj.getSubinterval(fieldName, level, 'subjects', subject, 'under', under);

                % Scenario (3)
            elseif isfield(obj.dataset{subject}.(string(level)), under) && ...
                    isfield(obj.dataset{subject}.(string(level)).(under), fieldName)
                feature = obj.dataset{subject}.(string(level)).(under).(fieldName);

                % Scenario (4)
            else
                obj = obj.(functionName)('subjects', subject, 'levels', level, 'modalities', functionFlag);
                feature = obj.dataset{subject}.(string(level)).(under).(fieldName);
            end
        
        end
        
    end

end

