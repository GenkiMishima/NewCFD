// MicroAVS Application: 7.00
// Date: 25 7ŒŽ 2003 11:54:54
// DataNames: 1 density_fld
// BitField: 0
// Methods: 1 10
// Actives: 2	 FLDBounds2D FLDIsoline2D
// Comments:  
// End
MicroAVS.MicroAVS MicroAVS {
   MicroAVS.MAviewer3D MAviewer {
      viewer_params {
         top_props {
            material = {0.3,0.7,0.,
50.};
            inherit = 0;
         };
         renderer = 1;
         top_xform {
            ocenter = {0.134005,0.749585,
0.};
            dcenter = {0.134005,0.749585,
0.};
            center = {0.134005,0.749585,
0.};
            mat = {
               2.44797,0.,0.,0.,0.,2.44797,0.,0.,0.,0.,2.44797,0.,0.,0.,
0.,1.
            };
            xlate = {-0.134004,-0.749585,
0.};
         };
         norm_scale = 0.8000000119;
         render_width = 1024;
         render_height = 624;
         auto_norm = 0;
         visible = 0;
      };
      objs_in => {
         MAviewerEditors.label_objs,<-.FLDBounds2D.obj_out,
         <-.FLDIsoline2D.obj_out};
      MAviewerEditors {
         LabelEditor {
            Params {
               label_editor_params {
                  Label1 {
                     red = 0.;
                     green = 0.;
                     blue = 0.;
                     height = 25;
                     label = "Particle Temperature (K) , dp=5e-7";
                     x_pos = -0.13;
                     y_pos = 0.5;
                     style = 1;
                  };
               };
            };
         };
         TransformEditor {
            MAGDxform_sub_obj {
               MAGDxform_save_mdl {
                  &Xform;
               };
               MAGDxform_update_mdl {
                  &Xform;
               };
            };
            Params {
               transform_editor_params {
                  center = 1;
                  scale = 0.;
               };
            };
         };
         BackcolorEditor {
            set_init_params {
               active = 0;
            };
            get_init_params {
               active = 1;
            };
            Params {
               backcolor_editor_params {
                  red = 1.;
                  green = 1.;
                  blue = 1.;
               };
            };
         };
      };
      MAviewer3Editors {
         PropsEditor {
            Params {
               props_editor_params {
                  ambient = 0.3;
                  diffuse = 0.7;
                  specular = 0.400000006;
                  gloss = 12.;
                  transparency = 1.;
                  inherit = 0;
                  red = 0.;
                  green = 0.;
                  blue = 0.;
               };
            };
            set_init_params {
               active = 0;
            };
            get_init_params {
               active = 1;
            };
         };
         CameraEditor {
            set_init_params {
               active = 0;
            };
            get_init_params {
               active = 1;
            };
         };
         AxisGraphEditor {
            Params {
               axis_editor_params {
                  axis1 {
                     start = -1.5;
                     end = 1.768009067;
                     step = 0.6536018133;
                  };
                  axis2 {
                     start = 0.;
                     end = 1.49916923;
                     step = 0.2998338461;
                  };
                  axis3 {
                     start = 0.;
                     end = 0.;
                     step = 0.;
                  };
               };
            };
         };
      };
   };
   MicroAVS.MAfieldReader density_fld {
      Params {
         field_params {
            each_values = {
               {
                  max_value=7.48528862,min_value=0.0275633838,,,,
               }};
            all_min_value = 3276.665771;
            all_max_value = 3497.396973;
            crop_dwn_set = 1;
         };
         end = 1;
         set_step_range = 0;
      };
      filename = "C:\\Documents and Settings\\toru\\My Documents\\Multi-Phase Fluid Dynamics of SRM\\Numerical Analysis\\multi_axi2d_L\\density.fld";
   };
   MicroAVS.field_networks.FLDBounds2D FLDBounds2D<NEvisible=1> {
      fld_in => <-.density_fld.field;
      Params {
         params {
            visible = 0;
         };
      };
      DataObject {
         Props {
            material = {0.3,0.7,0.4,
12.};
            inherit = 0;
            col = {0.,0.,0.};
         };
      };
   };
   MicroAVS.field_networks.FLDIsoline2D FLDIsoline2D<NEvisible=1> {
      Params {
         params {
            visible = 0;
         };
         colmap_params {
            curmin = 3241.03125;
            curmax = 3497.350098;
            min_col_val = 3241.03125;
            max_col_val = 3497.350098;
         };
         isoline_params {
            linenum = 50;
            level_min = 3241.03125;
            level_max = 3497.35;
            style = 1;
            height = 14;
            decimal = 0;
            level_state = 1;
            red = 0.;
            green = 0.;
            ContourFrequency = 10;
         };
         downsize_params {
            downsize = 1;
            max = 13;
         };
      };
      UI {
         MAisolineUI {
            text_font {
               font_panel {
                  x = 132;
                  y = 174;
               };
            };
         };
      };
      fld_in => <-.density_fld.field;
   };
};
MicroAVS.MAviewer.viewer_params.auto_norm = 0;// Auto Normalize Off
