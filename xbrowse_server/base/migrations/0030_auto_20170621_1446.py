# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2017-06-21 14:46
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('base', '0029_auto_20170608_0404'),
    ]

    operations = [
        migrations.AlterField(
            model_name='family',
            name='short_description',
            field=models.TextField(blank=True, default=b''),
        ),
    ]
