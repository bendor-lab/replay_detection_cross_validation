function varargout = data_view(varargin)
% DATA_VIEW M-file for data_view.fig
%      DATA_VIEW, by itself, creates a new DATA_VIEW or raises the existing
%      singleton*.
%
%      H = DATA_VIEW returns the handle to a new DATA_VIEW or the handle to
%      the existing singleton*.
%
%      DATA_VIEW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATA_VIEW.M with the given input arguments.
%
%      DATA_VIEW('Property','Value',...) creates a new DATA_VIEW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before data_view_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to data_view_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help data_view

% Last Modified by GUIDE v2.5 21-Jan-2019 19:20:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @data_view_OpeningFcn, ...
                   'gui_OutputFcn',  @data_view_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before data_view is made visible.
function data_view_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to data_view (see VARARGIN)

% Choose default command line output for data_view
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes data_view wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = data_view_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function cluster_pathway_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_pathway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cluster_pathway as text
%        str2double(get(hObject,'String')) returns contents of cluster_pathway as a double


% --- Executes during object creation, after setting all properties.
function cluster_pathway_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cluster_pathway (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in x_axis_variable.
function x_axis_variable_Callback(hObject, eventdata, handles)
% hObject    handle to x_axis_variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns x_axis_variable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from x_axis_variable


% --- Executes on selection change in y_axis_variable.
function y_axis_variable_Callback(hObject, eventdata, handles)
% hObject    handle to y_axis_variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns y_axis_variable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from y_axis_variable


% --- Executes during object creation, after setting all properties.
function y_axis_variable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_axis_variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in cluster_list.
function cluster_list_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cluster_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cluster_list


% --- Executes during object creation, after setting all properties.
function cluster_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cluster_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in color_variable.
function color_variable_Callback(hObject, eventdata, handles)
% hObject    handle to color_variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns color_variable contents as cell array
%        contents{get(hObject,'Value')} returns selected item from color_variable


% --- Executes during object creation, after setting all properties.
function color_variable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to color_variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in all_time.
function all_time_Callback(hObject, eventdata, handles)
% hObject    handle to all_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function scroll_index_Callback(hObject, eventdata, handles)
% hObject    handle to scroll_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function scroll_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scroll_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in hold.
function hold_Callback(hObject, eventdata, handles)
% hObject    handle to hold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function min_time_Callback(hObject, eventdata, handles)
% hObject    handle to min_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_time as text
%        str2double(get(hObject,'String')) returns contents of min_time as a double


% --- Executes during object creation, after setting all properties.
function min_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_time_Callback(hObject, eventdata, handles)
% hObject    handle to max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_time as text
%        str2double(get(hObject,'String')) returns contents of max_time as a double


% --- Executes during object creation, after setting all properties.
function max_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in shift_stretch.
function shift_stretch_Callback(hObject, eventdata, handles)
% hObject    handle to shift_stretch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in all_clusters.
function all_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to all_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of all_clusters


% --- Executes on slider movement.
function zoom_index_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function zoom_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoom_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in sort.
function sort_Callback(hObject, eventdata, handles)
% hObject    handle to sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in play.
function play_Callback(hObject, eventdata, handles)
% hObject    handle to play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in hold_axis.
function hold_axis_Callback(hObject, eventdata, handles)
% hObject    handle to hold_axis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold_axis


% --- Executes on button press in x_pos.
function x_pos_Callback(hObject, eventdata, handles)
% hObject    handle to x_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of x_pos


% --- Executes on button press in y_pos.
function y_pos_Callback(hObject, eventdata, handles)
% hObject    handle to y_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of y_pos


% --- Executes on selection change in eeg.
function eeg_Callback(hObject, eventdata, handles)
% hObject    handle to eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns eeg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eeg


% --- Executes during object creation, after setting all properties.
function eeg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in velocity.
function velocity_Callback(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of velocity


% --- Executes on button press in direction.
function direction_Callback(hObject, eventdata, handles)
% hObject    handle to direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of direction


% --- Executes on button press in lin_pos.
function lin_pos_Callback(hObject, eventdata, handles)
% hObject    handle to lin_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lin_pos


% --- Executes on button press in est_pos.
function est_pos_Callback(hObject, eventdata, handles)
% hObject    handle to est_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of est_pos


% --- Executes on button press in eeg2.
function eeg2_Callback(hObject, eventdata, handles)
% hObject    handle to eeg2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eeg2



function session_Callback(hObject, eventdata, handles)
% hObject    handle to session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of session as text
%        str2double(get(hObject,'String')) returns contents of session as a double


% --- Executes during object creation, after setting all properties.
function session_CreateFcn(hObject, eventdata, handles)
% hObject    handle to session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in position2.
function position2_Callback(hObject, eventdata, handles)
% hObject    handle to position2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of position2



function delay_value_Callback(hObject, eventdata, handles)
% hObject    handle to delay_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delay_value as text
%        str2double(get(hObject,'String')) returns contents of delay_value as a double


% --- Executes during object creation, after setting all properties.
function delay_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delay_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in replay_events.
function replay_events_Callback(hObject, eventdata, handles)
% hObject    handle to replay_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of replay_events


% --- Executes on button press in filter_option.
function filter_option_Callback(hObject, eventdata, handles)
% hObject    handle to filter_option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filter_option


% --- Executes on button press in analysis_type.
function analysis_type_Callback(hObject, eventdata, handles)
% hObject    handle to analysis_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analysis_type


% --- Executes on button press in analysis1.
function analysis1_Callback(hObject, eventdata, handles)
% hObject    handle to analysis1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analysis1


% --- Executes on button press in analysis2.
function analysis2_Callback(hObject, eventdata, handles)
% hObject    handle to analysis2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analysis2


% --- Executes on button press in analysis3.
function analysis3_Callback(hObject, eventdata, handles)
% hObject    handle to analysis3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analysis3


% --- Executes on button press in circle.
function circle_Callback(hObject, eventdata, handles)
% hObject    handle to circle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of circle


% --- Executes on button press in inner.
function inner_Callback(hObject, eventdata, handles)
% hObject    handle to inner (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of inner


% --- Executes on button press in ignore_direction.
function ignore_direction_Callback(hObject, eventdata, handles)
% hObject    handle to ignore_direction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ignore_direction



function highlight_Callback(hObject, eventdata, handles)
% hObject    handle to highlight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highlight as text
%        str2double(get(hObject,'String')) returns contents of highlight as a double


% --- Executes during object creation, after setting all properties.
function highlight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highlight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function highlight2_Callback(hObject, eventdata, handles)
% hObject    handle to highlight2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of highlight2 as text
%        str2double(get(hObject,'String')) returns contents of highlight2 as a double


% --- Executes during object creation, after setting all properties.
function highlight2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to highlight2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function folder_Callback(hObject, eventdata, handles)
% hObject    handle to folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of folder as text
%        str2double(get(hObject,'String')) returns contents of folder as a double


% --- Executes during object creation, after setting all properties.
function folder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in multiunit_activity.
function multiunit_activity_Callback(hObject, eventdata, handles)
% hObject    handle to multiunit_activity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multiunit_activity


% --- Executes on button press in next_replay.
function next_replay_Callback(hObject, eventdata, handles)
% hObject    handle to next_replay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in previous_replay.
function previous_replay_Callback(hObject, eventdata, handles)
% hObject    handle to previous_replay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in stimulus1.
function stimulus1_Callback(hObject, eventdata, handles)
% hObject    handle to stimulus1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of stimulus1



function subplot_choice_Callback(hObject, eventdata, handles)
% hObject    handle to subplot_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subplot_choice as text
%        str2double(get(hObject,'String')) returns contents of subplot_choice as a double


% --- Executes during object creation, after setting all properties.
function subplot_choice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subplot_choice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stimulus2.
function stimulus2_Callback(hObject, eventdata, handles)
% hObject    handle to stimulus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns stimulus2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stimulus2


% --- Executes during object creation, after setting all properties.
function stimulus2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stimulus2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reconstruct_bin_width_Callback(hObject, eventdata, handles)
% hObject    handle to reconstruct_bin_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reconstruct_bin_width as text
%        str2double(get(hObject,'String')) returns contents of reconstruct_bin_width as a double


% --- Executes during object creation, after setting all properties.
function reconstruct_bin_width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reconstruct_bin_width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_time_Callback(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_time as text
%        str2double(get(hObject,'String')) returns contents of start_time as a double


% --- Executes during object creation, after setting all properties.
function start_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stop_time_Callback(hObject, eventdata, handles)
% hObject    handle to stop_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stop_time as text
%        str2double(get(hObject,'String')) returns contents of stop_time as a double


% --- Executes during object creation, after setting all properties.
function stop_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stop_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function plotname_Callback(hObject, eventdata, handles)
% hObject    handle to plotname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of plotname as text
%        str2double(get(hObject,'String')) returns contents of plotname as a double


% --- Executes during object creation, after setting all properties.
function plotname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in number_of_active_clusters.
function number_of_active_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to number_of_active_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of number_of_active_clusters


% --- Executes on button press in ripple_power.
function ripple_power_Callback(hObject, eventdata, handles)
% hObject    handle to ripple_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ripple_power


% --- Executes on button press in up_down_state.
function up_down_state_Callback(hObject, eventdata, handles)
% hObject    handle to up_down_state (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of up_down_state


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function play_speed_Callback(hObject, eventdata, handles)
% hObject    handle to play_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of play_speed as text
%        str2double(get(hObject,'String')) returns contents of play_speed as a double


% --- Executes during object creation, after setting all properties.
function play_speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to play_speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
