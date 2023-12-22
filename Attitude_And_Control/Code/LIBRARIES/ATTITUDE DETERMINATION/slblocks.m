function blkStruct = slblocks
% This function specifies that the library 'mylib'
% should appear in the Library Browser with the 
% name 'My Library'

    Browser.Library = 'STATIC_ATT_DET';
    % 'mylib' is the name of the library

    Browser.Name = 'Attitude Determination (static)';
    % 'My Library' is the library name that appears
    % in the Library Browser

    blkStruct.Browser = Browser;
