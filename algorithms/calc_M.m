function [ M ] = calc_M( m_1, m_2 )
%helping function to calculate the Matrix M for the Helmke algorithm.
%   INPUT: Matched points of frame1 and frame2.

H = zeros(length(m_1),9);

for i=1:length(m_1)
    h = (m_1(i,:)'*m_2(i,:));
    H(i,:) = h(:)';
end

M = (1/length(m_1))* (H' * H);

end

