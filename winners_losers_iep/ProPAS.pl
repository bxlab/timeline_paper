#!/Users/guritsk1/miniconda2/bin/perl

#2013.7.22 add protein half-life calculation using the N-terminal rule.

use strict;
#require TK;
#use Tk;  
#use Tk::NoteBook;    
#use Tk::widgets qw/Dialog/; 
#use Tk::FileSelect;  
#use Tk::DialogBox;
#use Tk::Balloon;    
#use Tk::clipboard; 
#require Tk::ROText; 
##require 
#my %open_file;
#my %save_file;
#my $main = MainWindow->new;	#打开一个主窗口
##
##$main->minsize(qw(1024 768));	#定义最小窗口大小
##$main->maxsize(qw(600 600));	#定义最大窗口大小
##
#$main->configure(-title=>"ProPAS version 1.1", -background=>'blue');	#设置题目的名称和整个面板的颜色。可以用的有：red，blue，cyan
#$main->geometry('+100+50');	#设置打开时主窗口在屏幕上的相对位置
#
##定义主菜单
#my $menu_bar=$main->Frame(-relief=>'groove', -borderwidth=>3, -background=>'purple')->pack(-side=>'top', -fill=>'x');
#my $file_mb=$menu_bar->Menubutton(-text=>'File', -background => 'purple', -activebackground => 'cyan', -foreground=>'white')->pack(-side=>'left');
#$file_mb->command(-label=>'Exit', -activebackground=>'magenta', -command=>sub{$main->destroy});	
#$file_mb->separator();
#my $help_mb=$menu_bar->Menubutton(-text=>'Help', -background => 'purple', -activebackground => 'cyan', -foreground=>'white')->pack(-side=>'left');
#$help_mb->command(-label=>'About', -activebackground=>'magenta', -command=>\&about_txt);	
#$help_mb->command(-label=>'Help', -activebackground=>'magenta', -command=>\&help_txt);	
##About
#$help_mb->separator();
#    my $ver_dialog =   $main->Dialog(-title   => 'About ProPAS V1.1', 
#                                -text    => "Protein property analysis software (ProPAS) Version 1.1\n\nCopyright \251 2011 Beijing Proteome Research Center \n(Bioinformatics Group)\nAll rights reserved.", 
#                                -buttons => ['OK']);
#                               # -bitmap  => 'info'); 
#    my $menu = $help_mb->cget('-menu'); 
#    $menu->entryconfigure('About', -command => [$ver_dialog   => 'Show']); #这about要和前面的-label的一致。
##Help
#$help_mb->separator();
#    my $help_dialog =   $main->Dialog(-title   => 'ProPAS Help', 
#                                -text    => "ProPAS is used to calculate the basic property of peptides/proteins. The properties include the Isoelectric point (pI), Mass weight (MW) and Hydrophobicity (Hy). Additionally, amino acid composition of all peptides/proteins in the list was also calculated.\n\nDetail information could be found at http://bioinfor.hupo.org.cn/tools/ProPAS\n\nContact: songfengwu\100126.com\n", 
#                                -buttons => ['OK']);
#                               # -bitmap  => 'info'); 
##\n\nThere are three options for the input data: 1. Fasta format: it is the most popular format for protein sequences. Each protein is started with ‘>’ symbol, followed by the sequence name and the comment (separated by the space or table key). And then the sequence was showed in the following line. \n2. Single plain sequence: only contain one protein sequence.\n3. Peptide list: list multiple peptide sequences, one peptide one line, without any other symbols or sequence names and so on.\n\nFor the data input, you can paste the data to the text box directly, or load the file from “File open…” below the text box. If both of them above contain data/file, the pasted data would be ignored and the file would be used for the calculation.\n\nThe results could be exported to the file (“Save to file:” in the left bottom) or the right text box. If the output file directory and name were provided, the results would only be exported to the file but not the right text box.\n\nIn the “Parameter” section, you can determine which parameter would be calculated and statistic with special range and step.\n\nThe “Run” button is in the right bottom which will start the calculation. If the result is showed in the right text box, you can click “save to file…” in the bottom to save the results to the file, or click “Copy to clipboard” to save it in clipboard.\n
##$main->messageBox(-message => "Running complete. Total cost time: $costtime s.", -type => "ok"); #报错
#    my $menu2 = $help_mb->cget('-menu'); 
#    $menu2->entryconfigure('Help', -command => [$help_dialog   => 'Show']); #这about要和前面的-label的一致。
#
##开始
#my $book = $main->NoteBook()->pack( -fill=>'both', -expand=>1 );
##定义一个板
#my $tab1 = $book->add( "Sheet 1", -label=>"pI/MW/Hy/aaCom");
#my $table1 = $tab1->Frame(-relief=>'groove', -height=>50, -borderwidth=>3)->pack(-side=>'left', -fill=>'x');
#my $table1_1 = $table1->Frame(-relief=>'groove', -borderwidth=>3)->pack(-side=>'left', -fill=>'x');
##定义一个Frame
#my $vu_frm1=$table1->Frame(-relief=>'groove', -borderwidth=>3)->pack(-side=>'top', -fill=>'x');
##标记所需要的注释
#my $vu_lb1=$vu_frm1->Label(-text=>'Enter peptide/protein sequences with format: ', -foreground=>'black')->pack(-side=>'left');
##分3种格式：一个是fasta格式，另外一个是flat format，还一个是peptide-line format
######添加对格式的判断－－开始
##my $vu_frm1_1=$table1->Frame(-relief=>'groove', -borderwidth=>3)->pack(-side=>'top', -fill=>'x');
###标记所需要的注释
##my $vu_lb1_1=$vu_frm1_1->Label(-text=>'Select the file format: ', -foreground=>'black')->pack(-side=>'left');
#
#my $fileformat = "fastaf";
#
#my $radio_frame = $vu_frm1->Frame()->pack(-side => "top");
#$radio_frame->Label(-text=>" ")->pack(-side => "left");
#my $radio_blue = $radio_frame->Radiobutton(-text => "Fasta", -value => "fastaf",
#                                           -variable=> \$fileformat)->pack(-side => "left");
#my $radio_yellow = $radio_frame->Radiobutton(-text => "Single plain sequence", -value => "singlef",
#                                             -variable=> \$fileformat)->pack(-side => "left");
#my $radio_red = $radio_frame->Radiobutton(-text => "Peptide list", -value => "peptidef",
#                                          -variable=> \$fileformat)->pack(-side => "left");
#
#
#
######添加对格式的判断－－结束
##输入框
#my $vu_frm2=$table1->Frame(-relief=>'groove', -borderwidth=>3)->pack(-side=>'top', -fill=>'x');
#my $input_seq="333yutuytu";   
#my $vu_en1=$vu_frm2->Text(-width=>71, -background=>'white')->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#
##tie *INPUT, ref $vu_en1, vu_en1;
##文件输入
#my $vu_frm3=$table1->Frame(-relief=>'groove', -borderwidth=>5)->pack(-side=>'top', -fill=>'x');
##my $vu_en1=$vu_frm3->Text(-width=>10, -background=>'white')->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
##my $vu_en1=$vu_frm3->BrowseEntry(-width=>60, -background=>'white')->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $start_dir="";
#my $vu_en2=$vu_frm3->FileSelect(-directory=>'c:\\');
##my $dialog = $vu_frm3->Dialog(-text => 'Save File?', -bitmap => 'question', -title => 'Save File Dialog', -default_button => 'Yes', -buttons => [qw/Yes No Cancel/]);
##my $d = $vu_frm3->DialogBox(-title => "Title", -buttons => ["OK", "Cancel"]);
##my $inputonly = $vu_frm3->InputO(-width=>3);
#$vu_frm3->Label(-text=>'Open from file:  ', -foreground=>'black')->pack(-side=>'left');
#my $open_file_name_from_entry="";
#my $vu_en3=$vu_frm3->Entry(-width=>30, -background=>'white', -textvariable=>\$open_file_name_from_entry)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#$vu_frm3->Label(-text=>' ', -foreground=>'black')->pack(-side=>'left');
##my $file_name2="";
#my $go=$vu_frm3->Button(-text => 'File open...', -command =>\&pp_file_open)->pack(-side=>'left', -anchor=>'w');
#
#my $vu_frm4=$table1->Frame(-relief=>'groove', -borderwidth=>2)->pack(-side=>'top', -fill=>'x');
#my $vu_lb2=$vu_frm4->Label(-text=>'Parameter:', -foreground=>'black')->pack(-side=>'left');
#####################     pI    ###################
#my $vu_frm5=$table1->Frame(-relief=>'groove', -borderwidth=>2)->pack(-side=>'top', -fill=>'x');
#my $pI_control=1;
#$vu_frm5->Checkbutton(-text => "pI ", -state=>"normal", -variable=>\$pI_control)->pack(-side => 'left',
#                                              -expand => 1);
my $pI_stat=1;
#$vu_frm5->Checkbutton(-text => "Stat: ", -variable=>\$pI_stat)->pack(-side => 'left',
#                                              -expand => 1);
#my $vu_lb3=$vu_frm5->Label(-text=>'From', -foreground=>'black')->pack(-side=>'left');
my $pI_from=0;
#my $vu_en4=$vu_frm5->Entry (-width=>10, -background=>'white', -textvariable=>\$pI_from)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $vu_lb4=$vu_frm5->Label(-text=>'To', -foreground=>'black')->pack(-side=>'left');
#my $pI_to=14;
#my $vu_en5=$vu_frm5->Entry (-width=>10, -background=>'white', -textvariable=>\$pI_to)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $vu_lb5=$vu_frm5->Label(-text=>'Step', -foreground=>'black')->pack(-side=>'left');
#my $pI_step=0.05;
#my $vu_en6=$vu_frm5->Entry (-width=>10, -background=>'white', -textvariable=>\$pI_step)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#$vu_frm5->Label(-text=>'   ', -foreground=>'black')->pack(-side=>'left');
######################      MW         #########################
#my $vu_frm6=$table1->Frame(-relief=>'groove', -borderwidth=>2)->pack(-side=>'top', -fill=>'x');
#my $mw_control=1;
#$vu_frm6->Checkbutton(-text => "MW", -variable=>\$mw_control)->pack(-side => 'left',
#                                              -expand => 1);
my $mw_stat=1;
my $hl_stat=1;
#$vu_frm6->Checkbutton(-text => "Stat: ", -variable=>\$mw_stat)->pack(-side => 'left',
#                                              -expand => 1);
#my $vu_lb_mwfrom=$vu_frm6->Label(-text=>'From', -foreground=>'black')->pack(-side=>'left');
my $mw_from=0;
#my $vu_en_mwfrom=$vu_frm6->Entry(-width=>10, -background=>'white', -textvariable=>\$mw_from)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $vu_lb_mwto=$vu_frm6->Label(-text=>'To', -foreground=>'black')->pack(-side=>'left');
#my $mw_to=100000;
#my $vu_en_mwto=$vu_frm6->Entry(-width=>10, -background=>'white', -textvariable=>\$mw_to)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $vu_lb_mwstep=$vu_frm6->Label(-text=>'Step', -foreground=>'black')->pack(-side=>'left');
#my $mw_step=1000;
#my $vu_en_mwstep=$vu_frm6->Entry(-width=>10, -background=>'white', -textvariable=>\$mw_step)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#$vu_frm6->Label(-text=>'   ', -foreground=>'black')->pack(-side=>'left');
#
######################      Hy         #########################
#my $vu_frm7=$table1->Frame(-relief=>'groove', -borderwidth=>2)->pack(-side=>'top', -fill=>'x');
#my $hy_control=1;
#$vu_frm7->Checkbutton(-text => "Hy ", -variable=>\$hy_control)->pack(-side => 'left',
#                                              -expand => 1);
my $hy_stat=1;
#$vu_frm7->Checkbutton(-text => "Stat: ", -variable=>\$hy_stat)->pack(-side => 'left',
#                                              -expand => 1);
#my $vu_lb_hyfrom=$vu_frm7->Label(-text=>'From', -foreground=>'black')->pack(-side=>'left');
#my $hy_from=-2;
#my $vu_en_hyfrom=$vu_frm7->Entry	(-width=>10, -background=>'white', -textvariable=>\$hy_from)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $vu_lb_hyto=$vu_frm7->Label(-text=>'To', -foreground=>'black')->pack(-side=>'left');
#my $hy_to=1;
#my $vu_en_hyto=$vu_frm7->Entry	(-width=>10, -background=>'white', -textvariable=>\$hy_to)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#my $vu_lb_hystep=$vu_frm7->Label(-text=>'Step', -foreground=>'black')->pack(-side=>'left');
#my $hy_step=0.05;
#my $vu_en_hystep=$vu_frm7->Entry	(-width=>10, -background=>'white', -textvariable=>\$hy_step)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#$vu_frm7->Label(-text=>'   ', -foreground=>'black')->pack(-side=>'left');
#
######################      Amino acid composition         #########################
#
#my $vu_frm8=$table1->Frame(-relief=>'groove', -borderwidth=>2)->pack(-side=>'top', -fill=>'x');
my $aac_control=1;
#$vu_frm8->Checkbutton(-text => "Amino acid composition", -variable=>\$aac_control)->pack(-side => 'left',
#                                              -expand => 1);
#$vu_frm8->Label(-text=>'                                                                                                    ', -foreground=>'black')->pack(-side=>'left');
######################   计算和存储    #########################
#my $vu_frm9=$table1->Frame(-relief=>'groove', -borderwidth=>5)->pack(-side=>'top', -fill=>'x');
#$vu_frm9->Label(-text=>'Save to file:  ', -foreground=>'black')->pack(-side=>'left');
#my $save_file_name_from_entry="";
#my $vu_en7=$vu_frm9->Entry(-width=>30, -background=>'white', -textvariable=>\$save_file_name_from_entry)->pack(-side=>'left', -pady=>3);  #pady表示上下宽度
#$vu_frm9->Label(-text=>' ', -foreground=>'black')->pack(-side=>'left');
#my $go2=$vu_frm9->Button(-text => 'Browser', -command =>\&pp_file_save)->pack(-side=>'left', -anchor=>'w');
##注释
#my $vu_frm9_2=$table1->Frame(-relief=>'groove', -borderwidth=>2)->pack(-side=>'top', -fill=>'x');
#$vu_frm9_2->Label(-text=>"Note: if the file name with full directory is provided in the above box, the results would be only exported to \nthe file, and would not show in the right text box.                                                                                         ", -foreground=>'black')->pack(-side=>'left');
#
##my $run=$vu_frm9->Button(-text => 'Run' , -command =>\&calculate_phychem)->pack(-side=>'right', -anchor=>'w');
#
#
#
######################右边的输出框########################
#my $table1_2 = $tab1->Frame(-relief=>'groove', -height=>45, -borderwidth=>3)->pack(-side=>'left', -fill=>'x');
#my $out_frame=$table1_2->Frame(-relief=>'groove', -height=>42, -borderwidth=>3)->pack(-side=>'top', -fill=>'x');
#my $scrollbar = $out_frame->Scrollbar(-orient => 'vertical'); #->pack(-side=>'right');
#my $vu_en_out=$out_frame->ROText(-width=>70, -height=>42, -background=>'white', -yscrollcommand => ['set' , $scrollbar])->pack(-side=>'left', -pady=>5);  #pady表示上下宽度
#$scrollbar->configure(-command => ['yview', $vu_en_out]);
#$scrollbar->pack(-side => 'right', -fill => 'y');
#$table1_2->Label(-text=>'    ', -foreground=>'black')->pack(-side=>'left');
#my $run=$table1_2->Button(-text => 'Calculate' , -command =>\&calculate_phychem)->pack(-side=>'left', -anchor=>'w');
#$table1_2->Label(-text=>'            ', -foreground=>'black')->pack(-side=>'left');
#$table1_2->Label(-text=>'            ', -foreground=>'black')->pack(-side=>'right');
#my $go3=$table1_2->Button(-text => 'Copy to clipboard', -command =>\&pp_text_file_clipboard)->pack(-side=>'right', -anchor=>'w');
#$table1_2->Label(-text=>'            ', -foreground=>'black')->pack(-side=>'right');
#my $go4=$table1_2->Button(-text => 'Save to file...', -command =>\&pp_text_file_save)->pack(-side=>'right', -anchor=>'w');
#MainLoop;
#################图形化编码结束#####################


#print "$vu_en1\n";
#
#
##输入文件
#sub pp_file_open 
#{
##	my $types = [
##	    ['All Files',        '*',             ],
##        	];
##	my $filename ||= $vu_frm3->getOpenFile(-filetypes=>$types);
#	my $filename ||= $vu_frm3->getOpenFile();
#	#print "--------$filename---------\n";
#	$open_file{"phychem"}=$filename;
#	$vu_en3->delete(0.1, 'end');
#	$vu_en3->insert('end', "$filename");
#}
##输出文件
#sub pp_file_save 
#{
##	my $types = [
##	    ['All Files',        '*',             ],
##        	];
##	my $filename ||= $vu_frm9->getSaveFile(-filetypes=>$types);
#	my $filename ||= $vu_frm9->getSaveFile();
#	#print "--------$filename---------\n";
#	$save_file{"phychem"}=$filename;
#	$vu_en7->delete(0.1, 'end');
#	$vu_en7->insert('end', "$filename")
#	
#	
#}
#
#sub pp_text_file_save 
#{
##	my $types = [
##	    ['All Files',        '*',             ],
##        	];
##	my $filename ||= $table1_2->getSaveFile(-filetypes=>$types);
#	my $filename ||= $table1_2->getSaveFile();
#	#print "--------$filename---------\n";
#	#$save_file{"phychem"}=$filename;
#	$filename=~s/\\/\\\\/g if ($filename);
#	open (FILE_in, ">$filename")||die;
#	my $ROtext_output = $vu_en_out->get(0.1, 'end');
#	print FILE_in "$ROtext_output";
#	
#}
#
#sub  pp_text_file_clipboard
#{
#	#my $CLIP = Win32::Clipboard(); 
#	my $read_from_ROtext=$vu_en_out->get(0.1, 'end');
#	$go4->clipboardClear;
#	$go4->clipboardAppend("$read_from_ROtext");
#	#print $CLIP;
#	#$CLIP = $vu_en_out->get(0.1, 'end');
##	open (FILE_in, ">$filename")||die;
##	my $ROtext_output = $vu_en_out->get(0.1, 'end');
##	print FILE_in "$ROtext_output";
#	
#}

my $par;
#my $inputf;
my $fileformat;
my $open_file_name;
my $i;
my $save_file_name;
my $inputf;

my $pI_control=1;
#$pI_top=14;
my $pI_to=14;
my $pI_step=0.05;
my $mw_control=1;
my $hl_control=1;
my $mw_to=100000;
my $mw_step=1000;
my $hy_control=1;
my $hy_to=1;
my $hy_from=-2;
my $hy_step=0.05;
my $hl_method='m';
my $outfileformat='txt';
my $par_num=@ARGV;
#my $save_file_name2="";
#print "$par_num\n";
if ($par_num==0)
{
	print "ProPAS commandline version 1.02.\nCopyright @ Beijing Proteome Research Center.\n";
	print "More information could be found in the following website:\n";
	print 'http://bioinfo.hupo.org.cn/tools/ProPAS/propas.htm';
	print "\n";
	print 'Or contact songfengwu@126.com';
	print "\n\nIf you want to use the program, you should at lest add -i followed by the input file.\n\n";
	print "Usage:\n";
	print "Example: perl ProPAS.commandline.pl -f fasta -i input.fasta -o output.fasta\n\n";
	print "-if: format of the input file.\n  fasta(default):  Each protein is started with ‘>’ symbol, followed by\n    the sequence name and the comment (separated by the space or table key).\n    And then the sequence was showed in the following line.\n";
	print "  single: single plain sequence, only contain one protein sequence.\n";
	print "  peptide: peptide list, list multiple peptide sequences, one peptide one\n    line, without any other symbols or sequence names and so on.\n";
	print "-h: organisms for the half life calculation.\n";
	print "  m: Mammalian (default).\n";
	print "  y: Yeast.\n";
	print "  e: E. coli.\n";
	print "-i: input file name (Necessary).\n";
	print "-o: output file name. (Default: \$infile.propas.txt)\n";
	print "-of: format of ouput file.\n";
	print "  txt: the tab separated file format (default).\n";
	print "  csv: csv format that could be open by excel directly.\n";
#	die;

}
else
{
for ($i=0; $i<$par_num; $i++)
{
	$par=$ARGV[$i];
	if ($par eq '-if')
	{
		$i++;
		$par=$ARGV[$i];
		$fileformat=$par;
		if (($fileformat ne 'single') and ($fileformat ne 'peptide') and ($fileformat ne 'fasta'))
		{
			print "Error input for the input file format of \'-if\'\n";
			die;
		}
	}
	elsif ($par eq '-i')
	{
		$i++;
		$par=$ARGV[$i];
		$open_file_name=$par;
	}
	elsif ($par eq '-of')
	{
		$i++;
		$par=$ARGV[$i];
		$outfileformat=$par;
		if (($outfileformat ne 'txt') and ($outfileformat ne 'csv'))
		{
			print "Error input for the output file format of \'-of\'\n";
			die;
		}
	}
	elsif ($par eq '-h')
	{
		$i++;
		$par=$ARGV[$i];
		$hl_method=$par;
		if (($hl_method ne 'm') and ($hl_method ne 'y') and ($hl_method ne 'e'))
		{
			print "Error input for the organisms to calculate protein's half life\'-h\'\n";
			die;
		}
	}
	elsif ($par eq '-o')
	{
		$i++;
		$par=$ARGV[$i];
		$save_file_name=$par;
	}
	else
	{
		print "Not recognized parameter of \"$par\".\n";
		die;
	}
}
		

if (!($fileformat))  #  #fastaf，peptidef，singlef
{
	$fileformat='fasta';
}
if (!($open_file_name))  #  #fastaf，peptidef，singlef
{
	print "No input file, please enter the data for calculation.\n";
	die;
}
if (!($save_file_name))  #  #fastaf，peptidef，singlef
{
	$save_file_name="$open_file_name".'.propas.txt';
}

sub calculate_phychem;

calculate_phychem;

#$open_file_name="";
#$save_file_name="";
}

sub calculate_phychem
{
#######################	test
#	print "pI_control\t$widget->{pI_control}\n";
	#print "========$pI_control\n";
	#print "========pI_from\t$pI_from\n";
	#print "$input_seq\n";
	
	#print "$t444";
	
	
#######################	test
	my $starttime=time();
#	my $text_input = $vu_en1->get(0.1, 'end');
#	$text_input =~ s/\n(\n)+/\n/g;
#	my $text_input2=$text_input;
#	$text_input2 =~ s/[ \n]//g;
	
	#my $open_file_name=$open_file{"phychem"};
#	my $open_file_name=$open_file_name_from_entry;
	$open_file_name=~s/\\/\\\\/g if ($open_file_name);
	#my $save_file_name=$save_file{"phychem"};
#	my $save_file_name=$save_file_name_from_entry;
	$save_file_name=~s/\\/\\\\/g if ($save_file_name);
	if ($open_file_name)
	{
	#	open (DB_in, "$open_file_name")||die $main->messageBox(-message => "Could not open the file! \nMaybe the directory is not in english, or could not find the file!", -type => "ok");
		open (DB_in, "$open_file_name")||die;
	}
	else
	{
		#$main->messageBox(-message => "==$text_input2==!", -type => "ok");
		print "No inputfile\n";
		die;
#		if (!($text_input2))
#		{
#		 #	die $main->messageBox(-message => "Could not find input data!", -type => "ok");
#		 	print "Could not find input data!\n";
#		 	die;
#		}
	}
	if ($outfileformat eq 'txt')
	{
		
	}
	elsif ($outfileformat eq 'csv')
	{
	#	$save_file_name2=$save_file_name;
		if ($save_file_name=~/txt$/i)
		{
			$save_file_name=~s/txt$/csv/i;
		}
		else
		{
			$save_file_name.='.csv';
		}
	#	print "$save_file_name\n";
	}
	open (FileOUT, ">$save_file_name");
	#open (OUT_class_iep, ">protein_iep_class");
	#open (OUT_class_mw, ">protein_mw_class");
	#open (OUT_class_kd, ">protein_kd_class");
	#open (OUT_aa, ">aa_cont") || die ("Cannot open aa_cont");
	print "Running, plese waiting...";
#	$/ = "\n>"; # breaks at a \n> sequence (instead of \n)
	my $out_to_text="";
	if ($save_file_name)
	{
		if ($outfileformat eq 'txt')
		{
			print FileOUT "=====Result of the physicochemical property for each protein======\n\n" if (($hy_control==1) or ($mw_control==1) or ($hl_control==1) or ($pI_control==1));
			print FileOUT "Pro_id" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
			print FileOUT "\tHydrophobicity" if ($hy_control==1);
			print FileOUT "\tMW" if ($mw_control==1);
			print FileOUT "\tHalf_life" if ($hl_control==1);
			print FileOUT "\tpI" if ($pI_control==1);
			print FileOUT "\n" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
		}
		elsif ($outfileformat eq 'csv')
		{
			print FileOUT "\"=====Result of the physicochemical property for each protein======\"\n\n" if (($hy_control==1) or ($mw_control==1) or ($hl_control==1) or ($pI_control==1));
			print FileOUT "\"Pro_id\"" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
			print FileOUT ",\"Hydrophobicity\"" if ($hy_control==1);
			print FileOUT ",\"MW\"" if ($mw_control==1);
			print FileOUT ",\"Half_life\"" if ($hl_control==1);
			print FileOUT ",\"pI\"" if ($pI_control==1);
			print FileOUT "\n" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
		}
	}
	else
	{
		$out_to_text.="======Result of the physicochemical property for each protein======\n\n" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
		$out_to_text.="Pro_id" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
		$out_to_text.="\tHydrophobicity" if ($hy_control==1);
		$out_to_text.="\tMW" if ($mw_control==1);
		$out_to_text.="\tHalf_life" if ($hl_control==1);
		$out_to_text.="\tpI" if ($pI_control==1);
		$out_to_text.="\n" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
	}
	my $j=0;
	my $i=0;
	my @pi_count="";
	my @mw_count="";
	my @kd_count="";
	if ($pI_control==1)
	{
		my $pI_top=($pI_to)/$pI_step;
		for ($i=0; $i<=$pI_top; $i++)
		{
			$pi_count[$i]=0;
		}
	}
	$j=0;
	if ($mw_control==1)
	{
		my $mw_top=($mw_to)/$mw_step;
		for ($i=0; $i<=$mw_top; $i++)
		{
			$mw_count[$i]=0;
		}
	}
	if ($hy_control==1)
	{
		my $hy_top=($hy_to-$hy_from)/$hy_step;
		for ($i=0; $i<=$hy_top; $i++)
		{
			$kd_count[$i]=0;
		}
	}
	my $sum_aa = 0;
	my $sum_pi=0;
	my $sum_mw=0;
	my $sum_kd=0;
	my @sum_caa="";
	for ($i=0; $i<26; $i++)
	{
		$sum_caa[$i]=0;
	}
	
	my $file_in;
	my @text_input_array="";
	#if (!($text_input2))
	if (($open_file_name))
	{
		$file_in=<DB_in>;
	}
	else
	{
		print "No input file\n";
		die;
#		@text_input_array=split("\n", $text_input);
#		$file_in=$text_input_array[0];
#		$file_in=$file_in."\n";
	}
	#print "===File_in22\t$file_in\n====\n";
	
	my $i_text_array=0;
	my $name="";
	my $sequence="";
	while ($file_in) 
	{
		$name="";
		$sequence="";
		chomp $file_in;
		my $kd_protein;
		my $mw_protein;
		my $hl_protein;
		my $iep_protein;
		my $ii;
		my $jj;
		my $pI_top;
		my $mw_top;
		my $hy_top;
	#	if ($fileformat ne "fastaf")
	#	{
	#		print "Error: not fasta format--$fileformat\n";
	#	}
	#	print "test0: $file_in\n";
		#去除空行
		while (1)
		{
			chomp $file_in;
		#	print "$file_in\n";
			if ($file_in=~/[a-z\>]/i)
			{
		#		print "last\n";
				last;
			}
			else
			{
		#		print "next line: $file_in\n";
				if (($open_file_name))
				{
					$file_in=<DB_in>;
				}
				else
				{
					print "No input file\n";
					die;
				#	@text_input_array=split("\n", $text_input);
				#	$i_text_array++;
				#	$file_in=$text_input_array[$i_text_array];
				}
			}
		}
		
		
		if (($fileformat eq "fasta"))# and ($file_in=~/^>/))
		{
			if ($file_in!~/^>/)
			{
			#	$main->messageBox(-message => "Not Fasta format, re-input or reselect the file format please!", -type => "ok");
				print "Not Fasta format, re-input or reselect the file format please!\n";
				last;
			}
			#print "===File_in\t$file_in\n====\n";
			$name=$file_in;
			$name =~s/>//;
			$name =~s/ .*//;
			chomp $name;
			$sequence="";
			#print "===Name1\t$name\n====\n";
			if ($file_in=~/^>/)
			{
				if (($open_file_name))
				{
					$file_in=<DB_in>;
				}
				else
				{
					print "No input file\n";
					die;
#					$i_text_array++;
#					$file_in=$text_input_array[$i_text_array];
					#$file_in=$file_in."\n";
				}
			}
			while (($file_in !~ /^>/) and ($file_in))
			{
				chomp $file_in;
				$file_in=~s/ //g;
				$file_in=~s/\*//g;
			#	print "test: $file_in\n";
			#	print "$sequence==$file_in=\n\n";
				$sequence="$sequence"."$file_in" if ($file_in);
				if (($open_file_name))
				{
					$file_in=<DB_in>;
				}
				else
				{
					print "No input file\n";
					die;
					
#					$i_text_array++;
#					$file_in=$text_input_array[$i_text_array];
					#$file_in=$file_in."\n";
				}
			}
			
			$kd_protein = &kd_calc($sequence,1) if (($hy_control==1) or ($hy_stat==1));
			$mw_protein = &mw_calc($sequence,1) if (($mw_control==1) or ($mw_stat==1));
			$hl_protein = &half_life_calc($sequence,1) if (($hl_control==1) or ($hl_stat==1));
			$iep_protein = &iep_calc($sequence,1) if (($pI_control==1) or ($pI_stat==1));
			if ($save_file_name)  
			{
				if ($outfileformat eq 'txt')
				{
					#print "===Name\t$name\n====\n";
					print FileOUT "$name" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1));
					print FileOUT "\t$kd_protein" if ($hy_control==1);
					print FileOUT "\t$mw_protein" if ($mw_control==1);
					print FileOUT "\t$hl_protein" if ($hl_control==1);
					print FileOUT "\t$iep_protein" if ($pI_control==1);
					print FileOUT "\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1));
				}
				elsif ($outfileformat eq 'csv')
				{
					#print "===Name\t$name\n====\n";
					print FileOUT "\"$name\"" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
					print FileOUT ",\"$kd_protein\"" if ($hy_control==1);
					print FileOUT ",\"$mw_protein\"" if ($mw_control==1);
					print FileOUT ",\"$hl_protein\"" if ($hl_control==1);
					print FileOUT ",\"$iep_protein\"" if ($pI_control==1);
					print FileOUT "\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				}
			}
			else
			{
				#print "===Name\t$name\n====\n";
				$out_to_text.="$name" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1));
				$out_to_text.="\t$kd_protein" if ($hy_control==1);
				$out_to_text.="\t$mw_protein" if ($mw_control==1);
				$out_to_text.="\t$hl_protein" if ($hl_control==1);
				$out_to_text.="\t$iep_protein" if ($pI_control==1);
				$out_to_text.="\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
			}
			################ 统计等电点－start ####################
			if (($pI_stat==1))
			{
				chomp $iep_protein;
				$i=0;
				$sum_pi++;
				$pI_top=($pI_to)/$pI_step;
				
				for ($i=$pI_from; $i<=$pI_top; $i++)
				{
					if (($iep_protein>=($i*$pI_step)) and ($iep_protein<(($i+1)*$pI_step)))
					{
						$pi_count[$i]++;
					}
				}
			}
			################ 统计等电点－end   ####################
	        	
			################ 统计分子量－start ####################	
			if (($mw_stat==1))
			{
				chomp $mw_protein;
				#$i=0;
				$sum_mw++;
				$mw_top=($mw_to)/$mw_step;
				for ($i=$mw_from; $i<$mw_top; $i++)
				{
					$mw_count[$i]++ if (($mw_protein>=($i*$mw_step)) and ($mw_protein<(($i+1)*$mw_step)));
	        		
				}
				$mw_count[$mw_top]++ if ($mw_protein>=$mw_to);
			}
			################ 统计分子量－end   ####################
			
			################ 统计疏水性－start ####################	
			if (($hy_stat==1))
			{
				chomp $kd_protein;
				$i=0;
				$hy_top=($hy_to-$hy_from)/$hy_step;
				for ($i=1; $i<=$hy_top; $i++)
				{
					$kd_count[$i]++ if (($kd_protein>=(($i-1)*$hy_step+$hy_from)) and ($kd_protein<(($i)*$hy_step+$hy_from)));
					
				}
				$kd_count[$hy_top+1]++ if ($kd_protein>=$hy_to);
				$kd_count[0]++ if ($kd_protein<$hy_from);
				$sum_kd++;
			}
			################ 统计疏水性－end ######################	
	        	
			################ 统计氨基酸含量－start ####################	
			if ($aac_control==1)		
			{
				
				for ($i=65; $i<91; $i++)
				{
					$ii="";
					$jj="";
					$ii=chr("$i");
					$jj=chr($i+32);
					while ($sequence =~ /[$ii$jj]/g)
					{
						$sum_caa[$i-65]++; 
						$sum_aa++;
					}
				}
			}
			################ 统计氨基酸含量－end ######################	
		}
		elsif (($fileformat eq "peptide")) # and ($file_in=~/^>/))
		{
		#	print "test2: $file_in\n";
			#print "===File_in\t$file_in\n====\n";
		#	print "Peptide list format: $file_in\n";
			if ($file_in=~/^>/)
			{
				print "It might be Fasta format, reselect please!\n";
#				$main->messageBox(-message => "It might be Fasta format, reselect please!", -type => "ok"); #报错
				last;
			}
			
			chomp $file_in;
			$name=$file_in;
			$name =~s/>//;
			$name =~s/ .*//;
#			chomp $name;
			#print "===Name1\t$name\n====\n";
#			if (($open_file_name))
#			{
#				$file_in=<DB_in>;
#			}
#			else
#			{
#				$i_text_array++;
#				$file_in=$text_input_array[$i_text_array];
#				#$file_in=$file_in."\n";
#			}
			if (($file_in=~/[a-z]/i))# and ($file_in))
			{
				$sequence="$file_in";
				if (($open_file_name))
				{
					$file_in=<DB_in>;
				}
				else
				{
					print "No input file\n";
					die;
#					$i_text_array++;
#					$file_in=$text_input_array[$i_text_array];
					#$file_in=$file_in."\n";
				}
			}
			
			$kd_protein = &kd_calc($sequence,1) if (($hy_control==1) or ($hy_stat==1));
			$mw_protein = &mw_calc($sequence,1) if (($mw_control==1) or ($mw_stat==1));
			$hl_protein = &half_life_calc($sequence,1) if (($hl_control==1) or ($hl_stat==1));
			$iep_protein = &iep_calc($sequence,1) if (($pI_control==1) or ($pI_stat==1));
			if ($save_file_name)
			{
				if ($outfileformat eq 'txt')
				{
					#print "===Name\t$name\n====\n";
					print FileOUT "$name" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
					print FileOUT "\t$kd_protein" if ($hy_control==1);
					print FileOUT "\t$mw_protein" if ($mw_control==1);
					print FileOUT "\t$hl_protein" if ($hl_control==1);
					print FileOUT "\t$iep_protein" if ($pI_control==1);
					print FileOUT "\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				}
				elsif ($outfileformat eq 'csv')
				{
					#print "===Name\t$name\n====\n";
					print FileOUT "\"$name\"" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
					print FileOUT ",\"$kd_protein\"" if ($hy_control==1);
					print FileOUT ",\"$mw_protein\"" if ($mw_control==1);
					print FileOUT ",\"$hl_protein\"" if ($hl_control==1);
					print FileOUT ",\"$iep_protein\"" if ($pI_control==1);
					print FileOUT "\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				}
			}
			else
			{
				#print "===Name\t$name\n====\n";
				$out_to_text.="$name" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				$out_to_text.="\t$kd_protein" if ($hy_control==1);
				$out_to_text.="\t$mw_protein" if ($mw_control==1);
				$out_to_text.="\t$hl_protein" if ($hl_control==1);
				$out_to_text.="\t$iep_protein" if ($pI_control==1);
				$out_to_text.="\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
			}
			################ 统计等电点－start ####################
			if (($pI_stat==1))
			{
				chomp $iep_protein;
				$i=0;
				$sum_pi++;
				$pI_top=($pI_to)/$pI_step;
				
				for ($i=$pI_from; $i<=$pI_top; $i++)
				{
					$pi_count[$i]++ if (($iep_protein>=($i*$pI_step)) and ($iep_protein<(($i+1)*$pI_step)));
					
				}
			}
			################ 统计等电点－end   ####################
	        	
			################ 统计分子量－start ####################	
			if (($mw_stat==1))
			{
				chomp $mw_protein;
				#$i=0;
				$sum_mw++;
				$mw_top=($mw_to)/$mw_step;
				for ($i=$mw_from; $i<$mw_top; $i++)
				{
					$mw_count[$i]++ if (($mw_protein>=($i*$mw_step)) and ($mw_protein<(($i+1)*$mw_step)));
	        		
				}
				$mw_count[$mw_top]++ if ($mw_protein>=$mw_to);
			}
			################ 统计分子量－end   ####################
			
			################ 统计疏水性－start ####################	
			if (($hy_stat==1))
			{
				chomp $kd_protein;
				$i=0;
				$hy_top=($hy_to-$hy_from)/$hy_step;
				for ($i=1; $i<=$hy_top; $i++)
				{
					$kd_count[$i]++ if (($kd_protein>=(($i-1)*$hy_step+$hy_from)) and ($kd_protein<(($i)*$hy_step+$hy_from)));
					
				}
				$kd_count[$hy_top+1]++ if ($kd_protein>=$hy_to);
				$kd_count[0]++ if ($kd_protein<$hy_from);
				$sum_kd++;
			}
			################ 统计疏水性－end ######################	
	        	
			################ 统计氨基酸含量－start ####################	
			if ($aac_control==1)		
			{
				
				for ($i=0; $i<26; $i++)
				{
					$sum_caa[$i]=0;
				}
				for ($i=65; $i<91; $i++)
				{
					$ii="";
					$jj="";
					$ii=chr("$i");
					$jj=chr($i+32);
					while ($sequence =~ /[$ii$jj]/g)
					{
						$sum_caa[$i-65]++; 
						$sum_aa++;
					}
				}
			}
			################ 统计氨基酸含量－end ######################	
		}
		elsif (($fileformat eq "single"))# and ($file_in=~/^>/))
		{
		#	print "test3: $file_in\n";
			#print "===File_in\t$file_in\n====\n";
	#		$name=$file_in;
	#		$name =~s/>//;
	#		$name =~s/ .*//;
	#		chomp $name;
			$name="Single_seq";
			$sequence="";
			#print "===Name1\t$name\n====\n";
			if ($file_in=~/^>/)
			{
				print "It might be Fasta format, reselect please!\n";
#				$main->messageBox(-message => "It might be Fasta format, reselect please!", -type => "ok"); #报错
				last;
			}

	#		if ($file_in=~/^>/)
	#		{
	#			if (($open_file_name))
	#			{
	#				$file_in=<DB_in>;
	#			}
	#			else
	#			{
	#				$i_text_array++;
	#				$file_in=$text_input_array[$i_text_array];
	#				#$file_in=$file_in."\n";
	#			}
	#		}
			while (($file_in !~ /^>/) and ($file_in))
			{
				chomp $file_in;
				$sequence="$sequence"."$file_in";
				if (($open_file_name))
				{
					$file_in=<DB_in>;
				}
				else
				{
					print "No input file\n";
					die;
#					$i_text_array++;
#					$file_in=$text_input_array[$i_text_array];
					#$file_in=$file_in."\n";
				}
			}
			
			$kd_protein = &kd_calc($sequence,1) if (($hy_control==1) or ($hy_stat==1));
			$mw_protein = &mw_calc($sequence,1) if (($mw_control==1) or ($mw_stat==1));
			$hl_protein = &half_life_calc($sequence,1) if (($hl_control==1) or ($hl_stat==1));
			$iep_protein = &iep_calc($sequence,1) if (($pI_control==1) or ($pI_stat==1));
			if ($save_file_name)
			{
				if ($outfileformat eq 'txt')
				{
					#print "===Name\t$name\n====\n";
					print FileOUT "$name" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
					print FileOUT "\t$kd_protein" if ($hy_control==1);
					print FileOUT "\t$mw_protein" if ($mw_control==1);
					print FileOUT "\t$hl_protein" if ($hl_control==1);
					print FileOUT "\t$iep_protein" if ($pI_control==1);
					print FileOUT "\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				}
				elsif ($outfileformat eq 'csv')
				{
					#print "===Name\t$name\n====\n";
					print FileOUT "\"$name\"" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
					print FileOUT ",\"$kd_protein\"" if ($hy_control==1);
					print FileOUT ",\"$mw_protein\"" if ($mw_control==1);
					print FileOUT ",\"$hl_protein\"" if ($hl_control==1);
					print FileOUT ",\"$iep_protein\"" if ($pI_control==1);
					print FileOUT "\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				}
			}
			else
			{
				#print "===Name\t$name\n====\n";
				$out_to_text.="$name" if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
				$out_to_text.="\t$kd_protein" if ($hy_control==1);
				$out_to_text.="\t$mw_protein" if ($mw_control==1);
				$out_to_text.="\t$hl_protein" if ($hl_control==1);
				$out_to_text.="\t$iep_protein" if ($pI_control==1);
				$out_to_text.="\n"  if (($hy_control==1) or ($mw_control==1) or ($pI_control==1) or ($hl_control==1));
			}
			################ 统计等电点－start ####################
			if (($pI_stat==1))
			{
				chomp $iep_protein;
				$i=0;
				$sum_pi++;
				$pI_top=($pI_to)/$pI_step;
				
				for ($i=$pI_from; $i<=$pI_top; $i++)
				{
					$pi_count[$i]++ if (($iep_protein>=($i*$pI_step)) and ($iep_protein<(($i+1)*$pI_step)));
					
				}
			}
			################ 统计等电点－end   ####################
	        	
			################ 统计分子量－start ####################	
			if (($mw_stat==1))
			{
				chomp $mw_protein;
				#$i=0;
				$sum_mw++;
				$mw_top=($mw_to)/$mw_step;
				for ($i=$mw_from; $i<$mw_top; $i++)
				{
					$mw_count[$i]++ if (($mw_protein>=($i*$mw_step)) and ($mw_protein<(($i+1)*$mw_step)));
	        		
				}
				$mw_count[$mw_top]++ if ($mw_protein>=$mw_to);
			}
			################ 统计分子量－end   ####################
			
			################ 统计疏水性－start ####################	
			if (($hy_stat==1))
			{
				chomp $kd_protein;
				$i=0;
				$hy_top=($hy_to-$hy_from)/$hy_step;
				for ($i=1; $i<=$hy_top; $i++)
				{
					$kd_count[$i]++ if (($kd_protein>=(($i-1)*$hy_step+$hy_from)) and ($kd_protein<(($i)*$hy_step+$hy_from)));
					
				}
				$kd_count[$hy_top+1]++ if ($kd_protein>=$hy_to);
				$kd_count[0]++ if ($kd_protein<$hy_from);
				$sum_kd++;
			}
			################ 统计疏水性－end ######################	
	        	
			################ 统计氨基酸含量－start ####################	
			if ($aac_control==1)		
			{
				
				for ($i=0; $i<26; $i++)
				{
					$sum_caa[$i]=0;
				}
				for ($i=65; $i<91; $i++)
				{
					$ii="";
					$jj="";
					$ii=chr("$i");
					$jj=chr($i+32);
					while ($sequence =~ /[$ii$jj]/g)
					{
						$sum_caa[$i-65]++; 
						$sum_aa++;
					}
				}
			}
			################ 统计氨基酸含量－end ######################	
		}
		else
		{
		#	print "test4: $file_in\n";
			if (($open_file_name)) #跳过非“>”的
			{
				$file_in=<DB_in>;
			}
			else
			{
				print "No input file\n";
				die;
#				$i_text_array++;
	#			$file_in=$text_input_array[$i_text_array];
				#$file_in=$file_in."\n";
			}
		}
	}
	################ 统计结束，打印开始 ######################	

	######打印等电点统计数据#########
	my $iu;
	my $id;
	my $per=0;
	if ($pI_stat==1)
	{
		if ($save_file_name)
		{
			if  ($outfileformat eq 'txt')
			{
				print FileOUT "\n==========Result of the pI stat==========\n\n";
			}
		}
		else
		{
			$out_to_text.="\n==========Result of the pI stat==========\n\n";
		}
		my $pI_top=int (($pI_to)/$pI_step);
		my $pI_bottom=int ($pI_from/$pI_step);
		for ($i=$pI_bottom; $i<$pI_top; $i++)
		{
			
			$iu=($i+1)*$pI_step;
			$id=$i*$pI_step;
			$per=$pi_count[$i]/$sum_pi;
			$per*=10000;
			$per+=0.5;
			$per=~s/\..+//;
			$per/=10000;
			#print FileOUT "$id-$iu\t$pi_count[$i]\t$per\n";
			if ($save_file_name)
			{
				if ($outfileformat eq 'txt')
				{
					print FileOUT "$id-$iu\t$pi_count[$i]\t$per\n";
				}
				elsif ($outfileformat eq 'csv')
				{
					
				}
			}
			else
			{
				$out_to_text.="$id-$iu\t$pi_count[$i]\t$per\n";
			}
		}
	}
	#################################
	
	########打印分子量统计数据#######
	my $mwmw;
	if ($mw_stat==1)
	{
#
#for ($i=0; $i<50; $i++)
#{
#	$iu=2*($i+1)*1000;
#	$id=2*$i*1000;
#	$mwmw=$mw_count[2*$i]+$mw_count[2*$i+1];
#	##+$mw_count[2*$i+2]+$mw_count[2*$i+3];
#	$per=$mwmw/$sum_mw;
#	print OUT_class_mw "$id-$iu\t$mwmw\t$per\n";
#}
#$per=$mw_count[100]/$sum_mw;
#print OUT_class_mw " >100000\t$mw_count[100]\t$per\n";
#
		#print FileOUT "\n==========Result of the MW stat==========\n\n";
		if ($save_file_name)
		{
			if ($outfileformat eq 'txt')
			{
				print FileOUT "\n==========Result of the MW stat==========\n\n";
			}
		}
		else
		{
			$out_to_text.="\n==========Result of the MW stat==========\n\n";
		}
		my $mw_top=int (($mw_to)/$mw_step);
		my $mw_bottom=int ($mw_from/$mw_step);
		for ($i=$mw_bottom; $i<$mw_top; $i++)
		{
			$iu=($i+1)*$mw_step;
			$id=$i*$mw_step;
			$mwmw=$mw_count[$i];
			##+$mw_count[2*$i+2]+$mw_count[2*$i+3];
			$per=$mwmw/$sum_mw;
			$per*=10000;
			$per+=0.5;
			$per=~s/\..+//;
			$per/=10000;
			#print FileOUT "$id-$iu\t$mwmw\t$per\n";
			if ($save_file_name)
			{
				if ($outfileformat eq 'txt')
				{
					print FileOUT "$id-$iu\t$mwmw\t$per\n";
				}
				elsif ($outfileformat eq 'csv')
				{
					
				}
			}
			else
			{
				$out_to_text.="$id-$iu\t$mwmw\t$per\n";
			}
		}
		$per=$mw_count[$mw_top]/$sum_mw;
		$per*=10000;
		$per+=0.5;
		$per=~s/\..+//;
		$per/=10000;
		#print FileOUT " >$mw_to\t$mw_count[$mw_top]\t$per\n";
		if ($save_file_name)
		{
			if ($outfileformat eq 'txt')
			{
				print FileOUT " >$mw_to\t$mw_count[$mw_top]\t$per\n";
			}
			elsif  ($outfileformat eq 'csv')
			{
				
			}
		}
		else
		{
			$out_to_text.=" >$mw_to\t$mw_count[$mw_top]\t$per\n";
		}
	}
	#################################
	
	########打印疏水性统计数据######
	if ($hy_stat==1)
	{
		$per=$kd_count[0]/$sum_kd;
		$per*=10000;
		$per+=0.5;
		$per=~s/\..+//;
		$per/=10000;
		if ($save_file_name)
		{
			
			if  ($outfileformat eq 'txt')
			{
				print FileOUT "\n==========Result of the Hy stat==========\n\n";
				print FileOUT " <-2 \t$kd_count[0]\t$per\n";
			}
			elsif  ($outfileformat eq 'csv')
			{
				
			}
		}
		else
		{
			$out_to_text.="\n==========Result of the Hy stat==========\n\n";
			$out_to_text.=" <-2 \t$kd_count[0]\t$per\n";
		}			
		my $hy_top=int (($hy_to-$hy_from)/$hy_step);
		for ($i=1; $i<=$hy_top; $i++)
		{
			$iu=($i)*$hy_step+$hy_from;
			$id=($i-1)*$hy_step+$hy_from;
			$iu*=100;
			$iu+=0.5;
			$iu=~s/\..+//;
			$iu/=100;
			$id*=100;
			$id+=0.5;
			$id=~s/\..+//;
			$id/=100;
			$per=$kd_count[$i]/$sum_kd;
			$per*=10000;
			$per+=0.5;
			$per=~s/\..+//;
			$per/=10000;
			#print FileOUT "$id~$iu\t$kd_count[$i]\t$per\n";
			if ($save_file_name)
			{
				if ($outfileformat eq 'txt')
				{
					print FileOUT "$id~$iu\t$kd_count[$i]\t$per\n";
				}
				elsif ($outfileformat eq 'csv')
				{
					
				}
			}
			else
			{
				$out_to_text.="$id~$iu\t$kd_count[$i]\t$per\n";
			}
		}
		$per=$kd_count[$hy_top+1]/$sum_kd;
		$per*=10000;
		$per+=0.5;
		$per=~s/\..+//;
		$per/=10000;
		#print FileOUT " >1 \t$kd_count[$hy_top+1]\t$per\n";
		if ($save_file_name)
		{
			if ($outfileformat eq 'txt')
			{
				print FileOUT " >1 \t$kd_count[$hy_top+1]\t$per\n";
			}
			elsif ($outfileformat eq 'csv')
			{
				
			}
		}
		else
		{
			$out_to_text.=" >1 \t$kd_count[$hy_top+1]\t$per\n";
		}
	}
	#################################
	
	#######打印氨基酸含量统计数据###
	my $p_t=0;
	my $percent=0;
	my $kk;
	if ($aac_control==1)
	{
		#print FileOUT "\n==========Result of the amino acid composition stat==========\n\n";
		if ($save_file_name)
		{
			if ($outfileformat eq 'txt')
			{
				print FileOUT "\n==========Result of the amino acid composition stat==========\n\n";
			}
		}
		else
		{
			$out_to_text.="\n==========Result of the amino acid composition stat==========\n\n";
		}
		for ($i=0; $i<26; $i++)
		{
			$kk=chr($i+65);
			$percent=$sum_caa[$i]/$sum_aa if ($sum_aa!=0);
			$percent*=10000;
			$percent+=0.5;
			$percent=~s/\..+//;
			$percent/=10000;
			$p_t=$p_t+$percent;
#			print FileOUT "$kk\t$sum_caa[$i]\t$percent\n";
			if ($save_file_name)
			{
				if ($outfileformat eq 'txt')
				{
					print FileOUT "$kk\t$sum_caa[$i]\t$percent\n";
				}
				elsif ($outfileformat eq 'csv')
				{
					
				}
			}
			else
			{
				$out_to_text.="$kk\t$sum_caa[$i]\t$percent\n";
			}
		}
		#print FileOUT "Total sum = $sum_aa\n";
		if ($save_file_name)
		{
			if ($outfileformat eq 'txt')
			{
				print FileOUT "Total sum = $sum_aa\n";
			}
		}
		else
		{
			$out_to_text.="Total sum = $sum_aa\n";
		}
	}
	##################################
	 close(DB_in);
	 close FileOUT;
	 #close OUT_class_kd;
	 #close OUT_class_mw;
	 #close OUT_class_iep;
	 #close OUT_aa;
	print "\n================Finish!=======================";
	my $test_tmp="To be or not to be...\nThat is the question";
#	if (!($save_file_name))
#	{
	#	$vu_en_out->delete(0.1, 'end');	#清空文本输出框
	#	$vu_en_out->insert('end', "$out_to_text");	#从文本输出框输入
#	}	
#	my $endtime=time();
#	my $costtime=$endtime-$starttime;
#	print "Running complete. Total cost time: $costtime s.\n";
#	$main->messageBox(-message => "Running complete. Total cost time: $costtime s.", -type => "ok"); #报错
#	print "Total cost time: $costtime\n";
}
######################

 sub kd_calc {
    my $peptide    = shift;
    my $kd       = shift || 0;;
    my $hy_kd=0;
    my $index;
    my %aas;
if ($kd){ #Kyte-Doolittle
	%aas = (
A => 1.8,
R => -4.5,
N => -3.5,
D => -3.5,
C => 2.5,
Q => -3.5,
E => -3.5,
G => -0.4,
H => -3.2,
I => 4.5,
L => 3.8,
K => -3.9,
M => 1.9,
F => 2.8,
P => -1.6,
S => -0.8,
T => -0.7,
W => -0.9,
Y => -1.3,
V => 4.2
		 );
 }

 my $len = length($peptide);

    for ($index = 0; $index < $len; $index++) {

	my $letter = substr($peptide, $index, 1);
	$letter=uc($letter);

	$hy_kd += $aas{$letter};

    }
$hy_kd=$hy_kd/$len if ($len!=0);
$hy_kd*=10000;
$hy_kd+=0.5;
$hy_kd=~s/\..+//;
$hy_kd/=10000;

    return $hy_kd;    #疏水性计算
}
#######################

 sub mw_calc {
    my $peptide    = shift;
    my $kd       = shift || 0;;
     my $mass=0;
    my $index;
    my %aas;
if ($kd){ #Kyte-Doolittle
	%aas = (
A	=> 71.0788,
B => 114.6686,   # Asx   Aspartic acid or Asparagine 
C	=> 103.1388,
D	=> 115.0886,
E	=> 129.1155,
F	=> 147.1766,
G	=> 57.0519,
H	=> 137.1411,
I	=> 113.1594,
K	=> 128.1741,
L	=> 113.1594,
M	=> 131.1926,
N	=> 114.1038,
O	=> 237.3018,
P	=> 97.1167,
Q	=> 128.1307,
R	=> 156.1875,
S	=> 87.0782,
T	=> 101.1051,
U	=> 150.0388,
V	=> 99.1326,
W	=> 186.2132,
X => 111.1138, # Xaa   Any amino acid
Y	=> 163.176,
Z => 128.7531    #Glx   Glutamine or Glutamic acid
		 );
 }
 my $len = length($peptide);

    for ($index = 0; $index < $len; $index++) {

	my $letter = substr($peptide, $index, 1);
	$letter=uc($letter);

	$mass += $aas{$letter};

    }
$mass=$mass+18.015;
$mass+=0.0005;
$mass*=1000;
$mass=~s/\..+//;
$mass/=1000;
    return $mass;
}

#######################计算半衰期

 sub half_life_calc {
    my $peptide    = shift;
    my $kd       = shift || 0;;
     my $mass=0;
    my $index;
    my %aas;

#parameter for half life, depend on the N-aa
# Amino acid  Mammalian        Yeast     E. coli
# Ala  A        4.4 hour     >20 hour    >10 hour
# Arg  R          1 hour       2 min       2 min
# Asn  N        1.4 hour       3 min     >10 hour
# Asp  D        1.1 hour       3 min     >10 hour
# Cys  C        1.2 hour     >20 hour    >10 hour
# Gln  Q        0.8 hour      10 min     >10 hour
# Glu  E          1 hour      30 min     >10 hour
# Gly  G         30 hour     >20 hour    >10 hour
# His  H        3.5 hour      10 min     >10 hour
# Ile  I         20 hour      30 min     >10 hour
# Leu  L        5.5 hour       3 min       2 min
# Lys  K        1.3 hour       3 min       2 min
# Met  M         30 hour     >20 hour    >10 hour
# Phe  F        1.1 hour       3 min       2 min
# Pro  P        >20 hour     >20 hour      ?
# Ser  S        1.9 hour     >20 hour    >10 hour
# Thr  T        7.2 hour     >20 hour    >10 hour
# Trp  W        2.8 hour       3 min       2 min
# Tyr  Y        2.8 hour      10 min       2 min
# Val  V        100 hour     >20 hour    >10 hour
#if the first amino acid in the sequence is 'M', the amino acid would be removed.
    

if ($kd){ #Kyte-Doolittle
if ($hl_method eq 'm')
{
	%aas = (
A	=> 4.4,
C	=> 1.2,
D	=> 1.1,
E	=> 1,
F	=> 1.1,
G	=> 30,
H	=> 3.5,
I	=> 20,
K	=> 1.3,
L	=> 5.5,
M	=> 30,
N	=> 1.4,
P	=> 20,
Q	=> 0.8,
R	=> 1,
S	=> 1.9,
T	=> 7.2,
V	=> 100,
W	=> 2.8,
Y	=> 2.8
		 );
		}
elsif ($hl_method eq 'y')
{
	%aas = (
A	=> 20,
C	=> 20,
D	=> 3,
E	=> 30,
F	=> 3,
G	=> 20,
H	=> 10,
I	=> 30,
K	=> 3,
L	=> 3,
M	=> 20,
N	=> 3,
P	=> 20,
Q	=> 10,
R	=> 2,
S	=> 20,
T	=> 20,
V	=> 20,
W	=> 3,
Y	=> 10
		 );
	
}
elsif ($hl_method eq 'e')
{
	%aas = (
A	=> 10,
C	=> 10,
D	=> 10,
E	=> 10,
F	=> 2,
G	=> 10,
H	=> 10,
I	=> 10,
K	=> 2,
L	=> 2,
M	=> 10,
N	=> 10,
P	=> 0,
Q	=> 10,
R	=> 2,
S	=> 10,
T	=> 10,
V	=> 10,
W	=> 2,
Y	=> 2
	
		 );	
}
else
{
	print "Error: not recognize half life parameter\n";
	die;
}
 }
 my $len = length($peptide);

    for ($index = 0; $index < $len; $index++) {

	my $letter = substr($peptide, $index, 1);
	$letter=uc($letter);

	if (($letter ne 'M') or ($index>0))
	{ 
		$mass = $aas{$letter};
		last;
	}

    }
#$mass=$mass+18.015;
#$mass+=0.0005;
#$mass*=1000;
#$mass=~s/\..+//;
#$mass/=1000;
    return $mass;
}

###############################

sub iep_calc {


# Total number N of the amino acids Lys(K),Arg(R),His(H),Asp(D),Glu(E),Cys(C)and Tyr(Y) within the sequence.

my $N_K=0;
my $N_R=0;
my $N_H=0;
my $N_D=0;
my $N_E=0;
my $N_C=0;
my $N_Y=0;


my $sequence=shift;
while ($sequence =~ /K/ig){$N_K++}
while ($sequence =~ /R/ig){$N_R++}
while ($sequence =~ /H/ig){$N_H++}
while ($sequence =~ /D/ig){$N_D++}
while ($sequence =~ /E/ig){$N_E++}
while ($sequence =~ /C/ig){$N_C++}
while ($sequence =~ /Y/ig){$N_Y++}

# print fileout "K = $N_K\n";
# print fileout "R = $N_R\n";
# print fileout "H = $N_H\n";
# print fileout "D = $N_D\n";
# print fileout "E = $N_E\n";
# print fileout "C = $N_C\n";
# print fileout "Y = $N_Y\n\n";

# K values.

my $K_lys = 10**(-10.00);
my $K_arg = 10**(-12.00);
my $K_his = 10**(-5.98);
my $K_asp = 10**(-4.05);
my $K_glu = 10**(-4.45);
my $K_cys = 10**(-9.00);
my $K_tyr = 10**(-10.00);


my %pK_NT = (

A => 7.59,
C => 7.50,
D => 7.50,
E => 7.70,
F => 7.50,
G => 7.50,
H => 7.50,
I => 7.50,
K => 7.50,
L => 7.50,
M => 7.00,
N => 7.50,
P => 8.36,
Q => 7.50,
R => 7.50,
S => 6.93,
T => 6.82,
V => 7.44,
W => 7.50,
Y => 7.50
#B =>
#O =>
#U =>
#X =>
#Z =>
);

	my $NT_letter = substr ($sequence,0,1);
	$NT_letter=uc($NT_letter);
	my $K_NT = 0;
	$K_NT=10**(-$pK_NT{$NT_letter});

my %pK_CT = (

A => 3.55,
C => 3.55,
D => 4.55,
E => 4.75,
F => 3.55,
G => 3.55,
H => 3.55,
I => 3.55,
K => 3.55,
L => 3.55,
M => 3.55,
N => 3.55,
P => 3.55,
Q => 3.55,
R => 3.55,
S => 3.55,
T => 3.55,
V => 3.55,
W => 3.55,
Y => 3.55
);

	my $CT_letter = substr ($sequence,-1,1);
	$CT_letter=uc($CT_letter);
	my $K_CT = 10**(-$pK_CT{$CT_letter});
my $pH;
for ($pH=0;$pH<=14;$pH+=0.01){

	my $C_H = 10**(-$pH);

	my $CRp = $N_K*$C_H/($C_H+$K_lys) + $N_R*$C_H/($C_H+$K_arg) + $N_H*$C_H/($C_H+$K_his) + $C_H/($C_H+$K_NT);

	my $CRn = $N_D*$K_asp/($C_H+$K_asp) + $N_E*$K_glu/($C_H+$K_glu) + $N_C*$K_cys/($C_H+$K_cys) + $N_Y*$K_tyr/($C_H+$K_tyr) + $K_CT/($C_H+$K_CT);

	my $R = $CRp -$CRn;
	
	if (($pH==14) and ($R>=0))
	{
		return 14;
	}
	
	if ($R < 0){

	$pH -= 0.01;

$pH+=0.005;
$pH*=100;
$pH=~s/\..+//;
$pH/=100;

	return $pH;

		}
	}

}
