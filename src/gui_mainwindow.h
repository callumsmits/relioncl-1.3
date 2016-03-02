/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef GUI_MAINWINDOW_H_
#define GUI_MAINWINDOW_H_
#include "src/gui_jobwindow.h"
#include "src/gui_entries.h"

#define NR_BROWSE_TABS 13

#define DO_WRITE true
#define DONT_WRITE false
#define DO_READ true
#define DONT_READ false
#define DO_TOGGLE_CONT true
#define DONT_TOGGLE_CONT false
#define DO_GET_CL true
#define DONT_GET_CL false

// This class organises the main winfow of the relion GUI
static Fl_Hold_Browser *browser;
static Fl_Group        *browse_grp[NR_BROWSE_TABS];
static bool is_main_continue;
static GeneralJobWindow *job_general;
static CtffindJobWindow *job_ctffind;
static ManualpickJobWindow *job_manualpick;
static AutopickJobWindow *job_autopick;
static ExtractJobWindow *job_extract;
static SortJobWindow *job_sort;
static Class2DJobWindow *job_class2d;
static Class3DJobWindow *job_class3d;
static Auto3DJobWindow *job_auto3d;
static PostJobWindow *job_post;
static PolishJobWindow *job_polish;
static ResmapJobWindow *job_resmap;
static PublishJobWindow *job_publish;
//Toggle continue button
static Fl_Toggle_Button *toggle_continue;
// Run button
static Fl_Button *run_button;

class RelionMainWindow : public Fl_Window
{

public:

	// For Tabs
	Fl_Menu_Bar *menubar;
	Fl_Tabs *tabs;
	Fl_Group *tab0, *tab1, *tab2, *tab3, *tab4, *tab5;

    // Run button
    Fl_Button *print_CL_button, *cite_button, *display_button;

    // For job submission
    std::string outputname, final_command;
    std::vector<std::string> commands;

    //FileName for settings file
    std::string fn_settings;

	// Constructor with w x h size of the window and a title
	RelionMainWindow(int w, int h, const char* title);

    // Destructor
    ~RelionMainWindow(){};

    // Communicate with the different jobtype objects
    void jobCommunicate(bool do_write, bool do_read, bool do_toggle_continue, bool do_commandline, int this_job = 0);

private:


    // Vertical distance from the top
    int start_y;

    // Current height
    int current_y;


    /** Call-back functions for the Run button
     *  The method of using two functions of static void and inline void was copied from:
     *  http://www3.telus.net/public/robark/
     */

    static void cb_select_browsegroup(Fl_Widget*, void*);
    inline void cb_select_browsegroup_i();

    static void cb_toggle_continue(Fl_Widget*, void*);
    inline void cb_toggle_continue_i();

    static void cb_display(Fl_Widget*, void*);
    inline void cb_display_i();

    static void cb_run(Fl_Widget*, void*);
    inline void cb_run_i();

    static void cb_print_cl(Fl_Widget*, void*);
    inline void cb_print_cl_i();

    static void cb_menubar_load(Fl_Widget*, void*);
    inline void cb_menubar_load_i();

    static void cb_menubar_save(Fl_Widget*, void*);
    inline void cb_menubar_save_i();

    static void cb_menubar_reactivate_runbutton(Fl_Widget*, void*);
    inline void cb_menubar_reactivate_runbutton_i();

    static void cb_menubar_about(Fl_Widget*, void*);
    inline void cb_menubar_about_i();

    static void cb_menubar_quit(Fl_Widget*, void*);
    inline void cb_menubar_quit_i();
};

#endif /* GUI_MAINWINDOW_H_ */
