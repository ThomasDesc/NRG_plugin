from pymol import cmd


def output_message(output_box, text, type):
    # Type can be error, warning or valid
    out_color = None
    red = '<span style="color:red;">{}</span>'
    yellow = '<span style="color:orange;">{}</span>'
    green = '<span style="color:green;">{}</span>'
    if type == 'error':
        out_color = red
    elif type == 'warning':
        out_color = yellow
    elif type == ' valid':
        out_color = green
    output_box.append(out_color.format(text))


def refresh_dropdown(dropdown_to_refresh, output_box, filter_for='', no_warning=False):
    list_pymol_objects = cmd.get_names('all')
    if filter_for != '':
        list_pymol_objects = [x for x in list_pymol_objects if filter_for in x]
    if len(list_pymol_objects) == 0 and no_warning is False:
        output_message(output_box, 'No objects found', 'warning')
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(list_pymol_objects)
    if len(list_pymol_objects) > 0:
        dropdown_to_refresh.setCurrentText(list_pymol_objects[0])


def pymol_hide_structures(form):
    list_pymol_objects = cmd.get_names('all')
    list_pymol_objects = [x for x in list_pymol_objects if 'sph' in x]
    if form.button_hide.isChecked():
        form.button_hide.setText('Show')
        cmd.hide('everything', ','.join(list_pymol_objects))
    else:
        form.button_hide.setText('Hide')
        cmd.show('surface', ','.join(list_pymol_objects))


def get_mouse_config():
    config_mouse = ''
    try:
        name = cmd.get("button_mode_name")
        if name[0] == '1':
            config_mouse += 'one'
        elif name[0] == '2':
            config_mouse += 'two'
        elif name[0] == '3':
            config_mouse += 'three'
        config_mouse += '_button'
        if name[0] != '1':
            if name[9:] == 'Viewing':
                config_mouse += '_viewing'
            elif name[9:] == 'Editing':
                config_mouse += '_editing'
            elif name[9:] == 'Motions':
                config_mouse += '_motions'
        return config_mouse
    except:
        return 'three_button_viewing'

