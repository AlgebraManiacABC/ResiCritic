import sys
import tui
import files
import log
import state
from exit_codes import EXIT_SUCCESS, EXIT_FAILURE
from version import __version__

def main():

    tui.print_welcome()

    filename = files.get_excel_filename(sys.argv)
    if not filename:
        return EXIT_SUCCESS
    excel_df = files.parse_excel(filename)
    if excel_df is None or excel_df.empty:
        return EXIT_FAILURE

    log_filename = log.filename_to_log(filename)
    context = log.load_log(log_filename)
    if not context:
        return EXIT_FAILURE

    if context['step'] in state.state_functions or context['step'] == 'done':
        new_step = tui.already_analyzed_pick_step(filename,context)
        if new_step in state.state_functions:
            context['step'] = new_step
        else:
            context = log.init_log(log_filename)
            context['filename'] = filename
            context['dataframe'] = excel_df
            context['step'] = 'read-experiments'
    else:
        context = log.init_log(log_filename)
        context['filename'] = filename
        context['dataframe'] = excel_df
        context['step'] = 'read-experiments'

    while context['step'] in state.state_functions:
        err = state.run(context)
        if err != EXIT_SUCCESS:
            return EXIT_FAILURE
        log.save_log(context, log_filename)

    tui.done()
    return EXIT_SUCCESS

if __name__ == '__main__':
    sys.exit(main())
