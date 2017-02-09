function Intersections = cubeIntersect(Pt1,Pt2)
% Checks whether a line defined by two points intersects cubes.

Intersections = [];
% Initialise output.

Dvec = Pt2-Pt1;
% Distance between two points.

% Creating Bounding Box for Vessel
Imin = floor(min([Pt1(1),Pt2(1)]));
Imax = ceil(max([Pt1(1),Pt2(1)]));
Jmin = floor(min([Pt1(2),Pt2(2)]));
Jmax = ceil(max([Pt1(2),Pt2(2)]));
Kmin = floor(min([Pt1(3),Pt2(3)]));
Kmax = ceil(max([Pt1(3),Pt2(3)]));
% Limits the search of intersections to only those within the x, y
% and z limits of the two points.

for I2=Imin:Imax % Bounding X region of cylinder.
    for J2=Jmin:Jmax % Bounding Y region of cylinder.
        for K2=Kmin:Kmax % Bounding Z region of cylinder.
            
            boundbox_min = [I2-0.5,J2-0.5,K2-0.5];
            boundbox_max = [I2+0.5,J2+0.5,K2+0.5];
            % Values of cubes are assumed to be at centre so cube
            % walls are plus and minus 0.5 away.
            
            Intersect = false;
            % Initialise Intersect.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_min(1)-Pt1(1))/Dvec(1);
                if t2>0 && t2<1
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses lower x boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_min(2)-Pt1(2))/Dvec(2);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses lower y boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_min(3)-Pt1(3))/Dvec(3);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses lower z boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_max(1)-Pt1(1))/Dvec(1);
                if t2>0 && t2<1
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses upper x boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_max(2)-Pt1(2))/Dvec(2);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses upper y boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_max(3)-Pt1(3))/Dvec(3);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses upper z boundary.
            
            if Intersect
                Intersections(end+1,:) = [I2,J2,K2];
                % If the vessel intersected any of the boundries
                % then it intersects the voxel and it is added to
                % the list of intersections.
            end
        end
    end
end
end