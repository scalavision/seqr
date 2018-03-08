# -*- coding: utf-8 -*-
# Generated by Django 1.11 on 2018-03-07 18:03
from __future__ import unicode_literals

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('base', '0043_auto_20180219_1300'),
    ]

    operations = [
        migrations.CreateModel(
            name='AnalysedBy',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('date_saved', models.DateTimeField()),
                ('family', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='base.Family')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to=settings.AUTH_USER_MODEL)),
            ],
        ),
    ]