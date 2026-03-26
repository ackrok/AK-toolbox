function selection = overlayListDialog(targetFig, items, dlgTitle)
% overlayListDialog  Modal list dialog auto-sized to number and length of items
%   selection = overlayListDialog(targetFig, items, dlgTitle)
%   - targetFig : figure handle
%   - items     : cellstr or string array of options
%   - dlgTitle  : optional dialog title (default 'Select')
%   Returns selected item (char) or [] if cancelled.

    if nargin < 3, dlgTitle = 'Select'; end
    selection = [];

    % validate target figure
    if ~ishandle(targetFig) || ~strcmp(get(targetFig,'Type'),'figure')
        error('targetFig must be a valid figure handle.');
    end

    % normalize items to cellstr
    if isstring(items), items = cellstr(items); end

    % compute dialog size:
    % width: base + estimate from longest string (approx 7 px per char)
    maxChars = max(cellfun(@numel, items));
    W = max(220, 10 + round(maxChars * 7));
    W = min(W, 800); % clamp

    % height: header + listbox per item (20 px per item) + buttons area
    H = 40 + min(400, 20 * numel(items)) + 50;
    H = max(H, 120);

    % center dialog over target fig (pixel units)
    oldUnits = get(targetFig,'Units'); set(targetFig,'Units','pixels');
    pf = get(targetFig,'Position'); set(targetFig,'Units',oldUnits);
    left = round(pf(1) + (pf(3)-W)/2);
    bottom = round(pf(2) + (pf(4)-H)/2); % center vertically
    d = dialog('Name',dlgTitle,'Position',[left bottom W H]);

    % listbox (leave some margin)
    lbMargin = 10;
    lbHeight = H - 90;
    lb = uicontrol('Parent',d,'Style','listbox','String',items, ...
        'Position',[lbMargin 60 W-2*lbMargin lbHeight],'Max',1,'Min',0);

    % OK and Cancel buttons
    btnW = 80; btnH = 28;
    uicontrol('Parent',d,'Style','pushbutton','String','OK', ...
        'Position',[W-10-btnW 16 btnW btnH],'Callback',@okCB);
    uicontrol('Parent',d,'Style','pushbutton','String','Cancel', ...
        'Position',[10 16 btnW btnH],'Callback',@cancelCB);

    uiwait(d); % block until dialog closed

    function okCB(~,~)
        idx = lb.Value;
        if ~isempty(idx) && idx>=1
            s = lb.String;
            selection = char(s(idx));
        end
        delete(d);
    end

    function cancelCB(~,~)
        selection = [];
        delete(d);
    end
end
