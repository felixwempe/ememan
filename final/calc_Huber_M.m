function [ M_1, M_2, numb_good] = calc_Huber_M( E, m_1, m_2, delta )
%helping function to calculate the Matrix M for the Helmke algorithm.
%   INPUT: Matched points of frame1 and frame2. The Essential Matrix E and
%   a delta.
%   Output: Two Matrices where the data points are splitted. As well as
%   number of "good", i.e. m_1*E*m_2 < delta are "good" points otherwise 
%   the error is high.

H = zeros(length(m_1),9);

for i=1:length(m_1)
    h = (m_1(i,:)'*m_2(i,:));
    H(i,:) = h(:)';
end

Bool = abs(H * E(:)) <= delta;

if sum(Bool) == 0
    M_1 = zeros(9,9);
else
    M_1 = (1/sum(Bool))* (H(Bool,:)' * H(Bool,:));
end

if sum(~Bool == 0)
    M_2 = zeros(1,9);
else
    M_2 = H(~Bool,:);
end
    
numb_good = sum(Bool);

end

