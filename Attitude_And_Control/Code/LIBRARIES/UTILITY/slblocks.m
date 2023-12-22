function blkStruct = slblocks
% This function specifies that the library 'mylib'
% should appear in the Library Browser with the 
% name 'My Library'

    Browser.Library = 'my_math_ops';
    % 'mylib' is the name of the library

    Browser.Name = 'My Mathematics Operations';
    % 'My Library' is the library name that appears
    % in the Library Browser

    blkStruct.Browser = Browser;