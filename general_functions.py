from pymol import cmd


def refresh_dropdown(dropdown_to_refresh, filter_for=''):
    list_pymol_objects = cmd.get_names('all')
    if filter_for != '':
        list_pymol_objects = [x for x in list_pymol_objects if filter_for in x]
    dropdown_to_refresh.clear()
    dropdown_to_refresh.addItems(list_pymol_objects)
    if len(list_pymol_objects) > 0:
        dropdown_to_refresh.setCurrentText(list_pymol_objects[0])
