function [stpf] = mcstep(alphaLM, fLM, gLM, alphaTM, fTM, gTM, alphaUM, fUM, gUM)
    global isBracketed alphaMin alphaMax alphaU fU gU alphaL fL gL;
    isBound = false;
    derivativeSign = sign(gTM)*sign(gLM);

    % Check to see if Theorem 2.1 is violated
    if (gLM*(alphaTM-alphaLM)>=0.0)
        error("Line search failure: Interval cannot contain a minimizer!");        
    end
    
    % Case 1: Higher function value. Corresponds with case U1
    if (fTM>fLM)
		isBracketed = true;
		d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
		s = max([abs(d1) max([abs(gLM) abs(gTM)])]);
		d2 = s*sqrt((d1/s)*(d1/s) - (gLM/s)*(gTM/s));
		if (alphaTM<alphaLM) 
			d2 = -d2;
		end
		p = d1 + d2 - gLM;
		q = gTM + 2*d2 - gLM;
		cubic = alphaLM + p/q*(alphaTM-alphaLM);
		quadratic = alphaLM + (alphaTM-alphaLM)*(gLM/( ( (fLM-fTM)/(alphaTM-alphaLM) ) +gLM) )/2;
		if (abs(cubic-alphaLM)<abs(quadratic-alphaLM)) 
			stpf = cubic;
		else
			stpf = cubic + (quadratic-cubic)/2;
		end
    elseif (derivativeSign<0)
        % Case 2: lower function value with derivatives of opposite sign. 
        % Corresponds with case U3
		isBracketed = true;
		d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
		s = max([abs(d1) max([abs(gLM) abs(gTM)])]);
		d2 = s*sqrt((d1/s)*(d1/s) - (gLM/s)*(gTM/s));
		if (alphaTM>alphaLM) 
			d2 = -d2;
		end
		p = d1 + d2 - gTM;
		q = gLM + 2*d2 - gTM;
		cubic = alphaTM + p/q*(alphaLM - alphaTM);
		quadratic = alphaTM + (gTM/(gTM-gLM))*(alphaLM-alphaTM);
		if (abs(cubic-alphaTM)>abs(quadratic-alphaTM))
			stpf = cubic;
		else
			stpf = quadratic;
        end
    elseif abs(gTM)<abs(gLM)
        %Case 3: lower function value, derivatives of the same sign, with 
        %decreasing magnitude. Case U2
        isBound = true;
		d1 = gLM + gTM + 3*(fLM-fTM)/(alphaTM-alphaLM);
		s = max([abs(d1) max([abs(gLM) abs(gTM)])]);
        d2 = s*sqrt(max([0 ((d1/s)*(d1/s) - (gLM/s)*(gTM/s))]));
        if alphaTM>alphaLM
            d2 = -d2;
        end
        p = d1 + d2 - gTM;
		q = gLM + 2*d2 - gTM;
        if (p/q<0 && d2~=0)
            %If cubic tends to infinity in the direction of the step
			cubic = alphaTM + p/q*(alphaLM-alphaTM);
        elseif alphaTM>alphaLM
            % If the cubic estimation will be out of bounds, restrict it
			cubic = alphaMax;
        else
			cubic = alphaMin;
        end
		quadratic = alphaTM + (gTM/(gTM-gLM))*(alphaLM-alphaTM);
        if isBracketed
            if(abs(alphaTM-cubic) < abs(alphaTM-quadratic))
                stpf = cubic;
            else 
                stpf = quadratic;
            end
        else
            % Since the interval has not been bracketed, along with case U2 
            % the following conditions satisfy safeguarding conditions 
            % 2.2+2.3 and 2.4+2.5
            if (alphaTM>alphaLM)
                stpf = alphaMax;
            else
                stpf = alphaMin;
            end
        end
    else
        % Case 4: Lower function value, derivatives have same sign without 
        % decreasing magnitude. Case U2
        if isBracketed 
            d1 = gUM + gTM + 3*(fTM-fUM)/(alphaUM-alphaTM);
            s = max([abs(d1) max([abs(gUM) abs(gTM)])]);
            d2 = 2*sqrt((d1/s)*(d1/s) - (gUM/s)*(gTM/s));
            if (alphaTM>alphaUM)
                d2 = -d2;
            end
            p = d1 + d2 - gTM;
            q = gUM + 2*d2 - gTM;
            cubic = alphaTM + p/q*(alphaUM-alphaTM);
            stpf = cubic;
        elseif (alphaTM>alphaLM)
            % Since the interval has not been bracketed, along with case U2 
            % the following conditions satisfy safeguarding conditions 
            % 2.2+2.3 and 2.4+2.5
            stpf = alphaMax;
        else
			stpf = alphaMin;
        end
    end
    
    % Update interval of uncertainty (Updating algorithm)
    if fTM>fLM
        alphaU = alphaTM;
        fU = fTM;
		gU = gTM;
    else 
        if derivativeSign<0
            alphaU = alphaLM;
			fU = fLM;
			gU = gLM;
        end
		alphaL = alphaTM;
		fL = fTM;
		gL = gTM;
    end
    
    % Compute new safeguarded step
    stpf = min([max([stpf alphaMin]) alphaMax]);
    % Modified version of bounds in case3?
    if (isBracketed && isBound)	
        if alphaUM>alphaLM
            stpf = min([(alphaLM + 0.66*(alphaUM-alphaLM)) stpf]);
        else
            stpf = max([(alphaLM + 0.66*(alphaUM-alphaLM)) stpf]);
        end
    end
end